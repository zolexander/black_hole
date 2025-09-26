#include "engine.hpp"
#include <kerrintegrate.hpp>
#include <kerr_inline.hpp>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <optional>
#include "panel.hpp"
#include "blackholesim.hpp"
namespace BlackholeSim
{
    Engine::Engine() {}
    // Helperfunktionen
    std::optional<std::filesystem::path> Engine::abspath_no_traversal(
        const std::filesystem::path &basepath,
        const std::filesystem::path &relpath)
    {

        const auto abspath = std::filesystem::weakly_canonical(basepath / relpath);

        // thanks to https://stackoverflow.com/questions/1878001/how-do-i-check-if-a-c-stdstring-starts-with-a-certain-string-and-convert-a
        const auto index = abspath.string().rfind(basepath.string(), 0);
        if (index != 0)
        {
            return std::nullopt;
        }
        return abspath;
    }

    std::string Engine::readFromFile(const char *filePath)
    {
        if (!std::filesystem::exists(filePath))
        {
            std::cerr << "File does not exist: " << filePath << "\n";
            return "";
        }
        if (abspath_no_traversal(std::filesystem::current_path(), filePath) == std::nullopt)
        {
            std::cerr << "File path traversal detected: " << filePath << "\n";
            return "";
        }
        std::ifstream file;
        file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        std::string code;
        try
        {
            file.open(filePath);
            std::stringstream stream;
            stream << file.rdbuf();
            file.close();
            code = stream.str();
        }
        catch (std::ifstream::failure &e)
        {
            std::cerr << "ERROR::FILE_NOT_SUCCESSFULLY_READ: " << filePath << "\n";
        }
        return code;
    }

    bool Engine::initGL(const char *title = "Kerr photon demo", int w = 1200, int h = 800)
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
        if (!gladLoadGL())
        {
            std::cerr << "gladLoadGL failed\n";
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
        prog = Engine::loadShader("vs.frag", "fs.vert");
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

        // init photons
        resetPhotons(photonCount);

        return true;
    }

    void Engine::resetPhotons(int count)
    {
        photonCount = count;
        if (mode == Mode::Kerr)
        {
            kerrPhotons.clear();
            kerrPhotons.reserve(photonCount);
            // create L spread and start positions spread along vertical line at x0
            double x0 = -30.0;
            double yMin = -8.0, yMax = 8.0;
            for (int i = 0; i < photonCount; ++i)
            {
                double y0 = yMin + (yMax - yMin) * (i / double(std::max(1, photonCount - 1)));
                // convert to polar (r,phi)
                double r0 = std::sqrt(x0 * x0 + y0 * y0);
                double phi0 = std::atan2(y0, x0);
                double L = y0; // use y as impact parameter -> L
                Photon p(L);
                p.reset(L, {r0, phi0});
                p.trail.reserve(trailLength);
                kerrPhotons.push_back(std::move(p));
            }
        }
        else
        {
            testPhotons.clear();
            testPhotons.reserve(photonCount);
            double x0 = -30.0;
            double yMin = -8.0, yMax = 8.0;
            for (int i = 0; i < photonCount; ++i)
            {
                double y0 = yMin + (yMax - yMin) * (i / double(std::max(1, photonCount - 1)));
                TestPhoton p(x0, y0);
                p.reset(x0, y0);
                p.trail.reserve(trailLength);
                testPhotons.push_back(std::move(p));
            }
        }
    }
    // Update functions
    void Engine::updateKerrPhotons(double a_spin)
    {
        for (auto &p : kerrPhotons)
        {
            if (!p.alive)
                continue;
            // integrate in Kerr coords (toy integrator)
            p.s = integrateKerr(p.s, p.L, p.h, a_spin);
            // store cartesian point
            float x = (float)(p.s.r * cos(p.s.phi));
            float y = (float)(p.s.r * sin(p.s.phi));
            p.trail.emplace_back(x, y);
            // trim
            if (p.trail.size() > trailLength)
            {
                size_t remove = p.trail.size() - trailLength;
                p.trail.erase(p.trail.begin(), p.trail.begin() + remove);
            }
            // horizon test
            if (p.s.r < 2.0)
                p.alive = false;
        }
    }

    void Engine::updateTestPhotons(double h, double M)
    {
        for (auto &p : testPhotons)
        {
            if (!p.alive)
                continue;
            double r = std::sqrt(p.s.x * p.s.x + p.s.y * p.s.y);
            double rs = 2.0 * M;
            // gravitational acceleration (Newtonian central) as toy model
            double ax = -M * p.s.x / (r * r * r + 1e-9);
            double ay = -M * p.s.y / (r * r * r + 1e-9);
            p.s.vx += ax * h;
            p.s.vy += ay * h;
            double norm = std::sqrt(p.s.vx * p.s.vx + p.s.vy * p.s.vy);
            if (norm > 1e-12)
            {
                p.s.vx /= norm;
                p.s.vy /= norm;
            }
            p.s.x += p.s.vx * h * 10.0; // scale step for visibility
            p.s.y += p.s.vy * h * 10.0;
            p.trail.emplace_back((float)p.s.x, (float)p.s.y);
            if (p.trail.size() > trailLength)
            {
                size_t remove = p.trail.size() - trailLength;
                p.trail.erase(p.trail.begin(), p.trail.begin() + remove);
            }
            if (r < rs)
                p.alive = false;
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

            // update physics
            if (mode == Mode::Kerr)
            {
                updateKerrPhotons(a_spin);
            }
            else
            {
                updateTestPhotons(0.01, 1.0);
            }

            // debug output occasionally
            static int dbg = 0;
            if ((dbg++ % 300) == 0)
            {
                if (!kerrPhotons.empty())
                    std::cerr << "DBG[0] Kerr L=" << kerrPhotons.front().L
                              << " r=" << kerrPhotons.front().s.r
                              << " phi=" << kerrPhotons.front().s.phi
                              << " trail=" << kerrPhotons.front().trail.size()
                              << " alive=" << kerrPhotons.front().alive << std::endl;
                if (!testPhotons.empty())
                    std::cerr << "DBG[0] Test x=" << testPhotons.front().s.x
                              << " y=" << testPhotons.front().s.y
                              << " trail=" << testPhotons.front().trail.size()
                              << " alive=" << testPhotons.front().alive << std::endl;
            }

            // render
            int display_w, display_h;
            glfwGetFramebufferSize(window, &display_w, &display_h);
            glViewport(0, 0, display_w, display_h);
            glClearColor(0.02f, 0.02f, 0.03f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);

            // draw photons
            if (mode == Mode::Kerr)
                drawKerrPhotons(prog, proj, view);
            else
                drawTestPhotons(prog, proj, view);

            // ImGui render
            ImGui::Render();
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

            glfwSwapBuffers(window);
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