#ifndef BLACKHOLESIM_ENGINE_HPP
#define BLACKHOLESIM_ENGINE_HPP
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <random>
#include <chrono>
#include <sstream>

// OpenGL / GLFW / GLAD / ImGui / GLM includes
// Important: Include GLAD before any headers that might include OpenGL headers
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
    // Constants for better code readability
    namespace Constants {
        // Rendering constants
        constexpr float DEFAULT_ZOOM = 0.2f;
        constexpr double DEFAULT_SPIN = 0.7;
        constexpr int DEFAULT_PHOTON_COUNT = 200;
        constexpr size_t DEFAULT_TRAIL_LENGTH = 200;
        constexpr int TRAIL_SEGMENTS = 3;
        
        // Physics constants
        constexpr double SCHWARZSCHILD_RADIUS = 2.0;
        constexpr double INTEGRATION_STEP = 0.01;
        constexpr double MASS_DEFAULT = 1.0;
        constexpr double VELOCITY_SCALE = 10.0;
        constexpr double EPSILON = 1e-9;  // Small value to prevent division by zero
        constexpr double VELOCITY_THRESHOLD = 1e-12;
        
        // Initial conditions
        constexpr double INITIAL_X = -30.0;
        constexpr double Y_MIN = -8.0;
        constexpr double Y_MAX = 8.0;
        
        // Rendering colors
        constexpr float KERR_RED = 0.8f;
        constexpr float KERR_BLUE = 0.2f;
        constexpr float TEST_RED = 0.2f;
        constexpr float TEST_BLUE = 1.0f;
        constexpr float COLOR_GREEN = 0.2f;
        constexpr float TEST_COLOR_GREEN = 0.4f;
        
        // Alpha values for trail segments
        constexpr float ALPHA_FAINT = 0.05f;
        constexpr float ALPHA_MEDIUM = 0.25f;
        constexpr float ALPHA_BRIGHT = 1.0f;
        
        // Debug output frequency
        constexpr int DEBUG_OUTPUT_FREQUENCY = 300;
    }
    
    /**
     * @brief Configuration for photon trail rendering colors
     */
    struct PhotonRenderConfig {
        float red;
        float green;
        float blue;
    };
    class Engine
    {
        // Window & GL
    public:
        GLFWwindow *window = nullptr;

        // Photons
        std::vector<Photon> kerrPhotons;
        std::vector<TestPhoton> testPhotons;

        // Simulation parameters
        int photonCount = Constants::DEFAULT_PHOTON_COUNT;
        size_t trailLength = Constants::DEFAULT_TRAIL_LENGTH;
        double zoom = Constants::DEFAULT_ZOOM;
        double a_spin = Constants::DEFAULT_SPIN;
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
        // Physics update functions
        void updateKerrPhotons(double a_spin);
        void updateTestPhotons(double h, double M);
        
        // Rendering functions
        void drawKerrPhotons(GLuint shader, const glm::mat4 &projMat, const glm::mat4 &viewMat);
        void drawTestPhotons(GLuint shader, const glm::mat4 &projMat, const glm::mat4 &viewMat);
        
        // Photon initialization helpers
        void initializeKerrPhotons();
        void initializeTestPhotons();
        double calculatePhotonYPosition(int index, double yMin, double yMax) const;
        void trimPhotonTrail(std::vector<glm::vec2>& trail);
        
        // Generic rendering helper
        template<typename PhotonType>
        void renderPhotonTrails(GLuint shader, const glm::mat4& projMat, const glm::mat4& viewMat,
                               const std::vector<PhotonType>& photons, const PhotonRenderConfig& config);
        void renderTrailSegments(int vertexCount, const PhotonRenderConfig& config);
        void outputDebugInfo();
    };
}
#endif // BLACKHOLESIM_ENGINE_HPP