#ifndef BLACKHOLESIM_ENGINE_HPP
#define BLACKHOLESIM_ENGINE_HPP
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
#include <blackholesim.hpp>
#include <filesystem>
#include <optional>
namespace BlackholeSim
{
    class Engine
    {
        // Window & GL
    public:
        GLFWwindow *window = nullptr;

        // Photons
        std::vector<Photon> kerrPhotons;
        std::vector<TestPhoton> testPhotons;

        // Params
        int photonCount = 200;
        size_t trailLength = 200;
        double zoom = 0.2;
        double a_spin = 0.7;
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
        Engine();
        ~Engine();

        void run();
        bool initGL(const char *title, int w, int h);
        void resetPhotons(int count);

    private:
        std::optional<std::filesystem::path> abspath_no_traversal(
            const std::filesystem::path &basepath,
            const std::filesystem::path &relpath);
        std::string readFromFile(const char *filePath);
        GLuint loadShader(const char *vertexPath, const char *fragmentPath);
        // Konstruktor/Destruktor
        void updateKerrPhotons(double a_spin);
        void updateTestPhotons(double h, double M);
        void drawKerrPhotons(GLuint shader, const glm::mat4 &projMat, const glm::mat4 &viewMat);
        void drawTestPhotons(GLuint shader, const glm::mat4 &projMat, const glm::mat4 &viewMat);
    };
}
#endif // BLACKHOLESIM_ENGINE_HPP