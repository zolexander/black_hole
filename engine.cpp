#include "kerrintegrate.hpp"
#include "kerr_inline.hpp"
#include "engine.hpp"
#include <GLFW/glfw3.h>
#include <fstream>
#include <filesystem>
#include <optional>
#include "panel.hpp"
#include "blackhole_struct.hpp"
namespace BlackholeSim
{
    Engine::Engine() {}
    
    // ============================================================================
    // UTILITY FUNCTIONS
    // ============================================================================
    /**
     * @brief Validates that a relative path doesn't traverse outside the base directory
     * @param basepath The base directory path
     * @param relpath The relative path to validate
     * @return Optional absolute path if valid, nullopt if path traversal detected
     */
    std::optional<std::filesystem::path> Engine::abspath_no_traversal(const std::filesystem::path &basepath,
                                                                      const std::filesystem::path &relpath)
    {
        const auto abspath = std::filesystem::weakly_canonical(basepath / relpath);
        // Check if the resolved path starts with the base path to prevent directory traversal
        const auto index = abspath.string().rfind(basepath.string(), 0);
        if (index != 0)
        {
            return std::nullopt;
        }
        return abspath;
    }

    /**
     * @brief Safely reads a file's contents with path traversal protection
     * @param filePath Path to the file to read
     * @return File contents as string, empty string on error
     */
    std::string Engine::readFromFile(const char *filePath)
    {
        if (!filePath) {
            std::cerr << "ERROR: Null file path provided\n";
            return "";
        }
        
        if (!std::filesystem::exists(filePath)) {
            std::cerr << "ERROR: File does not exist: " << filePath << "\n";
            return "";
        }
        
        if (abspath_no_traversal(std::filesystem::current_path(), filePath) == std::nullopt) {
            std::cerr << "ERROR: File path traversal detected: " << filePath << "\n";
            return "";
        }
        
        try {
            std::ifstream file(filePath, std::ios::in | std::ios::binary);
            if (!file.is_open()) {
                std::cerr << "ERROR: Cannot open file: " << filePath << "\n";
                return "";
            }
            
            // Get file size for efficient memory allocation
            file.seekg(0, std::ios::end);
            const auto fileSize = file.tellg();
            file.seekg(0, std::ios::beg);
            
            std::string content;
            content.reserve(static_cast<size_t>(fileSize));
            content.assign(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>());
            
            return content;
        }
        catch (const std::exception& e) {
            std::cerr << "ERROR: Failed to read file " << filePath << ": " << e.what() << "\n";
            return "";
        }
    }

    // ============================================================================
    // INITIALIZATION FUNCTIONS
    // ============================================================================
    
    /**
     * @brief Initialize OpenGL context, shaders, and ImGui
     * @param title Window title
     * @param w Window width
     * @param h Window height
     * @return true if initialization successful, false otherwise
     */
    bool Engine::initGL(const char *title, int w, int h)
    {
        if (!glfwInit())
        {
            std::cerr << "GLFW init failed\n";
            return false;
        }

        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

        window = glfwCreateWindow(w, h, title, nullptr, nullptr);
        if (!window)
        {
            std::cerr << "glfwCreateWindow failed\n";
            glfwTerminate();
            return false;
        }

        glfwMakeContextCurrent(window);
        if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
        {
            std::cerr << "gladLoadGLLoader failed\n";
            return false;
        }

        // ImGui init
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImGuiIO &io = ImGui::GetIO();
        (void)io;
        ImGui::StyleColorsDark();
        ImGui_ImplGlfw_InitForOpenGL(window, true);
        ImGui_ImplOpenGL3_Init("#version 330");
        prog = Engine::loadShader("vs.vert", "fs.frag");
        colorLoc = glGetUniformLocation(prog, "color");
        zoomLoc = glGetUniformLocation(prog, "zoom");
        alphaLoc = glGetUniformLocation(prog, "alpha");
        projLoc = glGetUniformLocation(prog, "proj");
        viewLoc = glGetUniformLocation(prog, "view");

        // VAO/VBO setup
        glGenVertexArrays(1, &vao);
        glGenBuffers(1, &vbo);
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glEnableVertexAttribArray(0);
        // layout(location=0) = vec2
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (void *)0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);

        // blending
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE);
        glDisable(GL_DEPTH_TEST);

        // projection / view
        proj = glm::ortho(-40.0f, 40.0f, -30.0f, 30.0f, -1.0f, 1.0f);
        view = glm::mat4(1.0f);

        // Initialize photon simulation
        resetPhotons(photonCount);

        return true;
    }

    // ============================================================================
    // PHOTON MANAGEMENT
    // ============================================================================
    
    /**
     * @brief Reset and initialize photons for simulation
     * @param count Number of photons to create
     */
    void Engine::resetPhotons(int count)
    {
        photonCount = count;
        
        if (mode == Mode::Kerr) {
            initializeKerrPhotons();
        } else {
            initializeTestPhotons();
        }
    }
    
    /**
     * @brief Initialize Kerr photons with proper initial conditions
     */
    void Engine::initializeKerrPhotons()
    {
        kerrPhotons.clear();
        kerrPhotons.reserve(photonCount);
        
        // Create photons spread along vertical line at initial x position
        const double x0 = Constants::INITIAL_X;
        const double yMin = Constants::Y_MIN;
        const double yMax = Constants::Y_MAX;
        
        for (int i = 0; i < photonCount; ++i) {
            // Calculate y position with even distribution
            const double y0 = calculatePhotonYPosition(i, yMin, yMax);
            
            // Convert Cartesian to polar coordinates (r, phi)
            const double r0 = std::sqrt(x0 * x0 + y0 * y0);
            const double phi0 = std::atan2(y0, x0);
            const double L = y0; // Use y-coordinate as impact parameter L
            
            Photon photon(L);
            photon.reset(L, {r0, phi0});
            photon.trail.reserve(trailLength);
            kerrPhotons.emplace_back(std::move(photon));
        }
    }
    
    /**
     * @brief Initialize test photons for Newtonian simulation
     */
    void Engine::initializeTestPhotons()
    {
        testPhotons.clear();
        testPhotons.reserve(photonCount);
        
        const double x0 = Constants::INITIAL_X;
        const double yMin = Constants::Y_MIN;
        const double yMax = Constants::Y_MAX;
        
        for (int i = 0; i < photonCount; ++i) {
            const double y0 = calculatePhotonYPosition(i, yMin, yMax);
            
            TestPhoton photon(x0, y0);
            photon.reset(x0, y0);
            photon.trail.reserve(trailLength);
            testPhotons.emplace_back(std::move(photon));
        }
    }
    
    /**
     * @brief Calculate evenly distributed y-position for photon initialization
     * @param index Photon index
     * @param yMin Minimum y value
     * @param yMax Maximum y value
     * @return Calculated y position
     */
    double Engine::calculatePhotonYPosition(int index, double yMin, double yMax) const
    {
        const int denominator = std::max(1, photonCount - 1);
        return yMin + (yMax - yMin) * (index / static_cast<double>(denominator));
    }
    // ============================================================================
    // PHYSICS UPDATE FUNCTIONS
    // ============================================================================
    
    /**
     * @brief Update Kerr photons using relativistic integration
     * @param a_spin Black hole spin parameter
     */
    void Engine::updateKerrPhotons(double a_spin)
    {
        for (auto &photon : kerrPhotons) {
            if (!photon.alive) {
                continue;
            }
            
            // Integrate photon trajectory in Kerr coordinates
            photon.s = integrateKerr(photon.s, photon.L, photon.h, a_spin);
            
            // Convert polar coordinates to Cartesian for rendering
            const float x = static_cast<float>(photon.s.r * std::cos(photon.s.phi));
            const float y = static_cast<float>(photon.s.r * std::sin(photon.s.phi));
            photon.trail.emplace_back(x, y);
            
            // Maintain trail length by removing old points
            trimPhotonTrail(photon.trail);
            
            // Check if photon has crossed the event horizon
            if (photon.s.r < Constants::SCHWARZSCHILD_RADIUS) {
                photon.alive = false;
            }
        }
    }

    /**
     * @brief Update test photons using Newtonian gravity approximation
     * @param h Integration time step
     * @param M Mass of the central object
     */
    void Engine::updateTestPhotons(double h, double M)
    {
        for (auto &photon : testPhotons) {
            if (!photon.alive) {
                continue;
            }
            
            // Calculate distance from center
            const double r = std::sqrt(photon.s.x * photon.s.x + photon.s.y * photon.s.y);
            const double schwarzschildRadius = Constants::SCHWARZSCHILD_RADIUS * M;
            
            // Calculate Newtonian gravitational acceleration
            const double rCubed = r * r * r + Constants::EPSILON; // Prevent division by zero
            const double ax = -M * photon.s.x / rCubed;
            const double ay = -M * photon.s.y / rCubed;
            
            // Update velocity
            photon.s.vx += ax * h;
            photon.s.vy += ay * h;
            
            // Normalize velocity to light speed (c = 1 in our units)
            const double velocityMagnitude = std::sqrt(photon.s.vx * photon.s.vx + photon.s.vy * photon.s.vy);
            if (velocityMagnitude > Constants::VELOCITY_THRESHOLD) {
                photon.s.vx /= velocityMagnitude;
                photon.s.vy /= velocityMagnitude;
            }
            
            // Update position with scaled step for better visibility
            photon.s.x += photon.s.vx * h * Constants::VELOCITY_SCALE;
            photon.s.y += photon.s.vy * h * Constants::VELOCITY_SCALE;
            
            // Add point to trail
            photon.trail.emplace_back(static_cast<float>(photon.s.x), static_cast<float>(photon.s.y));
            
            // Maintain trail length
            trimPhotonTrail(photon.trail);
            
            // Check if photon has been absorbed
            if (r < schwarzschildRadius) {
                photon.alive = false;
            }
        }
    }
    
    /**
     * @brief Trim photon trail to maintain maximum length
     * @param trail Reference to the photon trail vector
     */
    void Engine::trimPhotonTrail(std::vector<glm::vec2>& trail)
    {
        if (trail.size() > trailLength) {
            const size_t removeCount = trail.size() - trailLength;
            trail.erase(trail.begin(), trail.begin() + removeCount);
        }
    }
    void Engine::drawKerrPhotons(GLuint shader, const glm::mat4 &projMat, const glm::mat4 &viewMat)
    {
        glUseProgram(shader);
        glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(projMat));
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(viewMat));
        glUniform1f(zoomLoc, (float)zoom);

        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);

        for (const auto &p : kerrPhotons)
        {
            if (p.trail.size() < 2)
                continue;
            // push trail data
            glBufferData(GL_ARRAY_BUFFER, p.trail.size() * sizeof(glm::vec2), p.trail.data(), GL_DYNAMIC_DRAW);

            // color based on frequency-shift-ish (we don't compute real freq shift here)
            float rcol = 0.8f;
            float bcol = 0.2f;

            int verts = (int)p.trail.size();
            const int segs = 3;
            for (int s = 0; s < segs; ++s)
            {
                int start = (verts * s) / segs;
                int end = (verts * (s + 1)) / segs;
                int count = end - start;
                if (count < 2)
                    continue;
                float alpha = (s == 0 ? 0.05f : (s == 1 ? 0.25f : 1.0f));
                glUniform3f(colorLoc, rcol, 0.2f, bcol);
                glUniform1f(alphaLoc, alpha);
                glDrawArrays(GL_LINE_STRIP, start, count);
            }
        }

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
        glUseProgram(0);
    }

    void Engine::drawTestPhotons(GLuint shader, const glm::mat4 &projMat, const glm::mat4 &viewMat)
    {
        // same code, but iterate testPhotons
        glUseProgram(shader);
        glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(projMat));
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(viewMat));
        glUniform1f(zoomLoc, (float)zoom);

        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);

        for (const auto &p : testPhotons)
        {
            if (p.trail.size() < 2)
                continue;
            glBufferData(GL_ARRAY_BUFFER, p.trail.size() * sizeof(glm::vec2), p.trail.data(), GL_DYNAMIC_DRAW);
            float rcol = 0.2f;
            float bcol = 1.0f;
            int verts = (int)p.trail.size();
            const int segs = 3;
            for (int s = 0; s < segs; ++s)
            {
                int start = (verts * s) / segs;
                int end = (verts * (s + 1)) / segs;
                int count = end - start;
                if (count < 2)
                    continue;
                float alpha = (s == 0 ? 0.05f : (s == 1 ? 0.25f : 1.0f));
                glUniform3f(colorLoc, rcol, 0.4f, bcol);
                glUniform1f(alphaLoc, alpha);
                glDrawArrays(GL_LINE_STRIP, start, count);
            }
        }

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
        glUseProgram(0);
    }

    GLuint Engine::loadShader(const char *vertexPath, const char *fragmentPath)
    {
        // 1. Read shader source
        std::string vertexCode = readFromFile(vertexPath);
        std::string fragmentCode = readFromFile(fragmentPath);
        const char *vShaderCode = vertexCode.c_str();
        const char *fShaderCode = fragmentCode.c_str();

        // 2. Compile vertex shader
        GLuint vertex = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertex, 1, &vShaderCode, NULL);
        glCompileShader(vertex);
        // Check for errors
        int success;
        char infoLog[512];
        glGetShaderiv(vertex, GL_COMPILE_STATUS, &success);
        if (!success)
        {
            glGetShaderInfoLog(vertex, 512, NULL, infoLog);
            std::cerr << "ERROR::VERTEX_SHADER::COMPILATION_FAILED\n"
                      << infoLog << "\n";
        }

        // 3. Compile fragment shader
        GLuint fragment = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragment, 1, &fShaderCode, NULL);
        glCompileShader(fragment);
        glGetShaderiv(fragment, GL_COMPILE_STATUS, &success);
        if (!success)
        {
            glGetShaderInfoLog(fragment, 512, NULL, infoLog);
            std::cerr << "ERROR::FRAGMENT_SHADER::COMPILATION_FAILED\n"
                      << infoLog << "\n";
        }

        // 4. Link shaders into a program
        GLuint program = glCreateProgram();
        glAttachShader(program, vertex);
        glAttachShader(program, fragment);
        glLinkProgram(program);
        glGetProgramiv(program, GL_LINK_STATUS, &success);
        if (!success)
        {
            glGetProgramInfoLog(program, 512, NULL, infoLog);
            std::cerr << "ERROR::PROGRAM::LINKING_FAILED\n"
                      << infoLog << "\n";
        }

        // 5. Cleanup shaders (already linked)
        glDeleteShader(vertex);
        glDeleteShader(fragment);

        return program;
    }

    void Engine::run()
    {
        if (!window)
            return;

        while (!glfwWindowShouldClose(window))
        {
            glfwPollEvents();
            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();
            Controlpanel(*this);

            // Update physics simulation
            if (mode == Mode::Kerr) {
                updateKerrPhotons(a_spin);
            } else {
                updateTestPhotons(Constants::INTEGRATION_STEP, Constants::MASS_DEFAULT);
            }

            // Periodic debug output for monitoring simulation state
            outputDebugInfo();

            // render
            int display_w, display_h;
            glfwGetFramebufferSize(window, &display_w, &display_h);
            glViewport(0, 0, display_w, display_h);
            glClearColor(0.02f, 0.02f, 0.03f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);
            // draw photons
            if (mode == Mode::Kerr) {
                drawKerrPhotons(prog, proj, view);
                if(showHorizons) {
                    double rp = bh.r_plus();
                    double rm = bh.r_minus();
                    if (std::isfinite(rp)) {
                       BlackholeSim::drawBlackHoleVisuals(bh, prog, vao, vbo, proj, view, (float)zoom);
                    }
                  
                }
            }
            // Finalize Dear ImGui for this frame (builds the command lists)
            ImGui::Render();

            // Reset GL state to what ImGui expects to avoid crashes in the backend
            glUseProgram(0);
            glBindVertexArray(0);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
            glDisable(GL_CULL_FACE);
            glDisable(GL_DEPTH_TEST);
            glEnable(GL_BLEND);
            glBlendEquation(GL_FUNC_ADD);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

            glfwSwapBuffers(window);
        }
    }

    /**
     * @brief Output periodic debug information about simulation state
     */
    void Engine::outputDebugInfo()
    {
        static int debugCounter = 0;
        if ((debugCounter++ % Constants::DEBUG_OUTPUT_FREQUENCY) == 0) {
            if (!kerrPhotons.empty()) {
                const auto& photon = kerrPhotons.front();
                std::cerr << "DEBUG[Kerr] L=" << photon.L
                          << " r=" << photon.s.r
                          << " phi=" << photon.s.phi
                          << " trail_size=" << photon.trail.size()
                          << " alive=" << (photon.alive ? "true" : "false") << std::endl;
            }
            if (!testPhotons.empty()) {
                const auto& photon = testPhotons.front();
                std::cerr << "DEBUG[Test] x=" << photon.s.x
                          << " y=" << photon.s.y
                          << " trail_size=" << photon.trail.size()
                          << " alive=" << (photon.alive ? "true" : "false") << std::endl;
            }
        }
    }

    Engine::~Engine()
    {
        if (vbo)
            glDeleteBuffers(1, &vbo);
        if (vao)
            glDeleteVertexArrays(1, &vao);
        if (prog)
            glDeleteProgram(prog);
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
        if (window)
        {
            glfwDestroyWindow(window);
            glfwTerminate();
        }
    }
}