// kerr_blackhole.cpp
// Minimal, self-contained example integrating ideas from the chat:
// - Two modes (Kerr / Test)
// - Photon trails as vector<glm::vec2>
// - Trail length controlled from ImGui
// - VAO/VBO rendering with shader (location = 0)
// - Additive blending glow + alpha per draw segment
//
// NOTE: Replace / integrate your real Kerr integrator (f_r/f_phi) if available.
// This demo uses a simplified "Kerr-like" integrator stub and a cartesian test integrator.

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <random>
#include <chrono>

// OpenGL / GLFW / GLAD / ImGui / GLM includes
#include <glad/glad.h>
#include <GLFW/glfw3.h>

// ImGui (assumes you have ImGui integrated)
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

// glm
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

////////////////////////////////////////////////////////////////////////////////
// Simple structures for Kerr and Cartesian states
////////////////////////////////////////////////////////////////////////////////

struct KerrState {
    double r;
    double phi;
};

struct CartesianState {
    double x, y;
    double vx, vy;
};

////////////////////////////////////////////////////////////////////////////////
// Photon class for Kerr-mode (stores KerrState, but trail is cartesian points)
////////////////////////////////////////////////////////////////////////////////

struct Photon {
    KerrState s;
    double L;
    double h;
    bool alive;
    std::vector<glm::vec2> trail;

    Photon(double L_ = 0.0) : L(L_), h(0.01), alive(true) {
        s = {4.0, M_PI};
        trail.reserve(256);
    }

    void reset(double L_, const KerrState &startState) {
        L = L_;
        s = startState;
        h = 0.01;
        alive = true;
        trail.clear();
        trail.reserve(256);
    }
};

////////////////////////////////////////////////////////////////////////////////
// TestPhoton class for cartesian-mode (keeps trail too)
////////////////////////////////////////////////////////////////////////////////

struct TestPhoton {
    CartesianState s;
    bool alive;
    std::vector<glm::vec2> trail;

    TestPhoton(double x0 = -20.0, double y0 = 0.0) {
        s = {x0, y0, 1.0, 0.0};
        alive = true;
        trail.reserve(256);
    }

    void reset(double x0, double y0) {
        s.x = x0; s.y = y0;
        s.vx = 1.0; s.vy = 0.0;
        alive = true;
        trail.clear();
        trail.reserve(256);
    }
};

////////////////////////////////////////////////////////////////////////////////
// Simple "integrator" for Kerr-like motion (stub).
// The real physics is more complicated — replace f_r/f_phi + integrator if available.
////////////////////////////////////////////////////////////////////////////////

inline double Delta(double r, double a) {
    return r*r - 2.0*r + a*a; // geometric units M=1 assumption (simplified)
}

// A sample R(r, L, a) placeholder; the real R_of_r is more complex.
inline double R_of_r(double r, double L, double a) {
    // Toy potential: allow radial motion; keep positive in some region.
    double U = 1.0 - 2.0/r;
    return std::max(0.0, U - (L*L)/(r*r));
}

// simplified dr/dlambda and dphi/dlambda functions (toy)
inline double f_r(double r, double L, double a) {
    double val = R_of_r(r, L, a);
    return -sqrt(std::max(0.0, val)) / (r*r + 1e-9);
}
inline double f_phi(double r, double L, double a) {
    double sig = r*r;
    double del = Delta(r, a);
    double A = (r*r + a*a) - a*L;
    // simplified variant that depends only on r,L,a
    return (a/del)*A - (a - L)/sig;
}

// RKF-like integrator for KerrState (toy RK4 with small step)
KerrState integrateKerr(const KerrState &s0, double L, double &h, double a_spin) {
    // simple RK4 integrator on (r,phi) with derivatives f_r, f_phi
    auto f = [&](const KerrState &st)->KerrState {
        return { f_r(st.r, L, a_spin), f_phi(st.r, L, a_spin) };
    };

    KerrState k1 = f(s0);
    KerrState s2  = { s0.r + 0.5*h*k1.r, s0.phi + 0.5*h*k1.phi };
    KerrState k2 = f(s2);
    KerrState s3  = { s0.r + 0.5*h*k2.r, s0.phi + 0.5*h*k2.phi };
    KerrState k3 = f(s3);
    KerrState s4  = { s0.r + h*k3.r, s0.phi + h*k3.phi };
    KerrState k4 = f(s4);

    KerrState out;
    out.r   = s0.r + (h/6.0)*(k1.r + 2*k2.r + 2*k3.r + k4.r);
    out.phi = s0.phi + (h/6.0)*(k1.phi + 2*k2.phi + 2*k3.phi + k4.phi);
    return out;
}

////////////////////////////////////////////////////////////////////////////////
// Engine: holds photons, testPhotons, mode, GUI parameters, VBO/VAO/shader
////////////////////////////////////////////////////////////////////////////////

enum class Mode { Kerr = 0, Test = 1 };

struct Engine {
    // Window & GL
    GLFWwindow* window = nullptr;

    // Photons
    std::vector<Photon> kerrPhotons;
    std::vector<TestPhoton> testPhotons;

    // Params
    int photonCount = 200;
    size_t trailLength = 200;
    double zoom = 0.2;
    double a_spin = 0.7; // initial spin for Kerr
    Mode mode = Mode::Kerr;

    // GL assets
    GLuint prog = 0;
    GLuint vao = 0;
    GLuint vbo = 0;
    GLint colorLoc = -1;
    GLint zoomLoc = -1;
    GLint alphaLoc = -1;
    GLint projLoc = -1;
    GLint viewLoc = -1;

    // camera projection
    glm::mat4 proj, view;

    // RNG for initial y distribution
    std::default_random_engine rng;
    std::uniform_real_distribution<double> dist{-5.0, 5.0};

    Engine() {}

    bool initGL(const char* title = "Kerr photon demo", int w = 1200, int h = 800) {
        if (!glfwInit()) {
            std::cerr << "GLFW init failed\n";
            return false;
        }

        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

        window = glfwCreateWindow(w, h, title, nullptr, nullptr);
        if (!window) {
            std::cerr << "glfwCreateWindow failed\n";
            glfwTerminate();
            return false;
        }

        glfwMakeContextCurrent(window);
        if (!gladLoadGL()) {
            std::cerr << "gladLoadGL failed\n";
            return false;
        }

        // ImGui init
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImGuiIO &io = ImGui::GetIO(); (void)io;
        ImGui::StyleColorsDark();
        ImGui_ImplGlfw_InitForOpenGL(window, true);
        ImGui_ImplOpenGL3_Init("#version 330");

        // simple shaders
        const char* vs = R"(
            #version 330 core
            layout(location = 0) in vec2 pos;
            uniform mat4 proj;
            uniform mat4 view;
            uniform float zoom;
            void main() {
                gl_Position = proj * view * vec4(pos * zoom, 0.0, 1.0);
            }
        )";
        const char* fs = R"(
            #version 330 core
            out vec4 FragColor;
            uniform vec3 color;
            uniform float alpha;
            void main() {
                FragColor = vec4(color, alpha);
            }
        )";

        prog = makeProgram(vs, fs);
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
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);
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

    void resetPhotons(int count) {
        photonCount = count;
        if (mode == Mode::Kerr) {
            kerrPhotons.clear();
            kerrPhotons.reserve(photonCount);
            // create L spread and start positions spread along vertical line at x0
            double x0 = -30.0;
            double yMin = -8.0, yMax = 8.0;
            for (int i = 0; i < photonCount; ++i) {
                double y0 = yMin + (yMax - yMin) * (i / double(std::max(1, photonCount - 1)));
                // convert to polar (r,phi)
                double r0 = std::sqrt(x0*x0 + y0*y0);
                double phi0 = std::atan2(y0, x0);
                double L = y0; // use y as impact parameter -> L
                Photon p(L);
                p.reset(L, {r0, phi0});
                p.trail.reserve(trailLength);
                kerrPhotons.push_back(std::move(p));
            }
        } else {
            testPhotons.clear();
            testPhotons.reserve(photonCount);
            double x0 = -30.0;
            double yMin = -8.0, yMax = 8.0;
            for (int i = 0; i < photonCount; ++i) {
                double y0 = yMin + (yMax - yMin) * (i / double(std::max(1, photonCount - 1)));
                TestPhoton p(x0, y0);
                p.reset(x0, y0);
                p.trail.reserve(trailLength);
                testPhotons.push_back(std::move(p));
            }
        }
    }

    // Update functions
    void updateKerrPhotons(double a_spin) {
        for (auto &p : kerrPhotons) {
            if (!p.alive) continue;
            // integrate in Kerr coords (toy integrator)
            p.s = integrateKerr(p.s, p.L, p.h, a_spin);
            // store cartesian point
            float x = (float)(p.s.r * cos(p.s.phi));
            float y = (float)(p.s.r * sin(p.s.phi));
            p.trail.emplace_back(x, y);
            // trim
            if (p.trail.size() > trailLength) {
                size_t remove = p.trail.size() - trailLength;
                p.trail.erase(p.trail.begin(), p.trail.begin() + remove);
            }
            // horizon test
            if (p.s.r < 2.0) p.alive = false;
        }
    }

    void updateTestPhotons(double h, double M) {
        for (auto &p : testPhotons) {
            if (!p.alive) continue;
            double r = std::sqrt(p.s.x*p.s.x + p.s.y*p.s.y);
            double rs = 2.0 * M;
            // gravitational acceleration (Newtonian central) as toy model
            double ax = -M * p.s.x / (r*r*r + 1e-9);
            double ay = -M * p.s.y / (r*r*r + 1e-9);
            p.s.vx += ax * h;
            p.s.vy += ay * h;
            double norm = std::sqrt(p.s.vx*p.s.vx + p.s.vy*p.s.vy);
            if (norm > 1e-12) {
                p.s.vx /= norm;
                p.s.vy /= norm;
            }
            p.s.x += p.s.vx * h * 10.0; // scale step for visibility
            p.s.y += p.s.vy * h * 10.0;
            p.trail.emplace_back((float)p.s.x, (float)p.s.y);
            if (p.trail.size() > trailLength) {
                size_t remove = p.trail.size() - trailLength;
                p.trail.erase(p.trail.begin(), p.trail.begin() + remove);
            }
            if (r < rs) p.alive = false;
        }
    }

    // Draw helpers
    void drawKerrPhotons(GLuint shader, const glm::mat4 &projMat, const glm::mat4 &viewMat) {
        glUseProgram(shader);
        glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(projMat));
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(viewMat));
        glUniform1f(zoomLoc, (float)zoom);

        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);

        for (const auto &p : kerrPhotons) {
            if (p.trail.size() < 2) continue;
            // push trail data
            glBufferData(GL_ARRAY_BUFFER, p.trail.size() * sizeof(glm::vec2), p.trail.data(), GL_DYNAMIC_DRAW);

            // color based on frequency-shift-ish (we don't compute real freq shift here)
            float rcol = 0.8f;
            float bcol = 0.2f;

            int verts = (int)p.trail.size();
            const int segs = 3;
            for (int s = 0; s < segs; ++s) {
                int start = (verts * s) / segs;
                int end = (verts * (s+1)) / segs;
                int count = end - start;
                if (count < 2) continue;
                float alpha = (s==0? 0.05f : (s==1? 0.25f : 1.0f));
                glUniform3f(colorLoc, rcol, 0.2f, bcol);
                glUniform1f(alphaLoc, alpha);
                glDrawArrays(GL_LINE_STRIP, start, count);
            }
        }

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
        glUseProgram(0);
    }

    void drawTestPhotons(GLuint shader, const glm::mat4 &projMat, const glm::mat4 &viewMat) {
        // same code, but iterate testPhotons
        glUseProgram(shader);
        glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(projMat));
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(viewMat));
        glUniform1f(zoomLoc, (float)zoom);

        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);

        for (const auto &p : testPhotons) {
            if (p.trail.size() < 2) continue;
            glBufferData(GL_ARRAY_BUFFER, p.trail.size() * sizeof(glm::vec2), p.trail.data(), GL_DYNAMIC_DRAW);
            float rcol = 0.2f;
            float bcol = 1.0f;
            int verts = (int)p.trail.size();
            const int segs = 3;
            for (int s = 0; s < segs; ++s) {
                int start = (verts * s) / segs;
                int end = (verts * (s+1)) / segs;
                int count = end - start;
                if (count < 2) continue;
                float alpha = (s==0? 0.05f : (s==1? 0.25f : 1.0f));
                glUniform3f(colorLoc, rcol, 0.4f, bcol);
                glUniform1f(alphaLoc, alpha);
                glDrawArrays(GL_LINE_STRIP, start, count);
            }
        }

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
        glUseProgram(0);
    }

    void run() {
        if (!window) return;

        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();

            // ImGui frame
            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();

            // Control window
            ImGui::Begin("Kontrolle");
            if (ImGui::SliderInt("Photonen", &photonCount, 10, 2000)) {
                resetPhotons(photonCount);
            }
            double amin = 0.0, amax = 1.0;
            if (ImGui::SliderScalar("Spin a/M", ImGuiDataType_Double, &a_spin, &amin, &amax)) {
                // nothing else needed here; integrator reads a_spin each step
            }
            size_t trmin = 10, trmax = 5000;
            uint64_t tmp = (uint64_t)trailLength;
            if (ImGui::SliderScalar("Trail", ImGuiDataType_U64, &tmp, &trmin, &trmax)) {
                trailLength = (size_t)tmp;
                // reserve for existing photons
                for (auto &p : kerrPhotons) p.trail.reserve(trailLength);
                for (auto &p : testPhotons) p.trail.reserve(trailLength);
            }
            ImGui::SliderFloat("Zoom", (float*)&zoom, 0.01f, 2.0f, "%.3f");
            if (ImGui::RadioButton("Kerr-Modus", (int*)&mode, (int)Mode::Kerr)) {
                resetPhotons(photonCount);
            }
            ImGui::SameLine();
            if (ImGui::RadioButton("Test-Modus", (int*)&mode, (int)Mode::Test)) {
                resetPhotons(photonCount);
            }

            if (ImGui::Button("Reset Photonen")) {
                resetPhotons(photonCount);
            }

            ImGui::Text("Photons: %d  TrailLen: %zu", photonCount, trailLength);
            ImGui::End();

            // update physics
            if (mode == Mode::Kerr) {
                updateKerrPhotons(a_spin);
            } else {
                updateTestPhotons(0.01, 1.0);
            }

            // debug output occasionally
            static int dbg = 0;
            if ((dbg++ % 300) == 0) {
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
            if (mode == Mode::Kerr) drawKerrPhotons(prog, proj, view);
            else drawTestPhotons(prog, proj, view);

            // ImGui render
            ImGui::Render();
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

            glfwSwapBuffers(window);
        }
    }

    // Shader helpers
    static GLuint compileShader(GLenum type, const char *src) {
        GLuint s = glCreateShader(type);
        glShaderSource(s, 1, &src, nullptr);
        glCompileShader(s);
        GLint ok = 0; glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
        if (!ok) {
            char buf[1024]; glGetShaderInfoLog(s, 1024, nullptr, buf);
            std::cerr << "Shader compile error: " << buf << std::endl;
        }
        return s;
    }

    static GLuint makeProgram(const char *vsSrc, const char *fsSrc) {
        GLuint vs = compileShader(GL_VERTEX_SHADER, vsSrc);
        GLuint fs = compileShader(GL_FRAGMENT_SHADER, fsSrc);
        GLuint p = glCreateProgram();
        glAttachShader(p, vs);
        glAttachShader(p, fs);
        glLinkProgram(p);
        GLint ok = 0; glGetProgramiv(p, GL_LINK_STATUS, &ok);
        if (!ok) {
            char buf[1024]; glGetProgramInfoLog(p, 1024, nullptr, buf);
            std::cerr << "Program link error: " << buf << std::endl;
        }
        glDeleteShader(vs);
        glDeleteShader(fs);
        return p;
    }

    ~Engine() {
        if (vbo) glDeleteBuffers(1, &vbo);
        if (vao) glDeleteVertexArrays(1, &vao);
        if (prog) glDeleteProgram(prog);
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
        if (window) {
            glfwDestroyWindow(window);
            glfwTerminate();
        }
    }
};

////////////////////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////////////////////

int main() {
    Engine engine;
    if (!engine.initGL("Photon Trails Demo", 1200, 800)) {
        return -1;
    }
    engine.run();
    return 0;
}
