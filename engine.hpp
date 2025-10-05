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
#include "blackhole_struct.hpp"
#include "utils.hpp"

namespace BlackholeSim
{
    // Constants for better code readability
    namespace Constants {
        // Rendering constants
        constexpr float DEFAULT_ZOOM = 0.9f;
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
        constexpr double INITIAL_X = -20.0;
        // Narrow band around photon critical impact parameter (~5.2 for M=1)
        constexpr double Y_MIN = -6.0;
        constexpr double Y_MAX = 6.0;
        
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
 

    struct Vertex {
        glm::vec2 pos;
        glm::vec3 color;
    };

    class Engine
    {
        // Window & GL
    public:
        GLFWwindow *window = nullptr;
        BlackHole bh = BlackHole();

        bool showHorizons = true;
        bool showErgosphere = true;
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
        // GPU photon path state (compute)
        bool useGPUPaths = false;
        BlackholeSim::Utils::Shader photonCompute; // compute_photons.comp
        GLuint ssboPhotons = 0; // binding = 0
        GLuint ssboTrails = 0;  // binding = 1
        int gpuTrailLen = (int)Constants::DEFAULT_TRAIL_LENGTH;
        unsigned int frameIndex = 0;
        Engine();
        ~Engine();
        Utils::Shader shader1;
        void run();
        bool initGL(const char *title, int w, int h);
        void resetPhotons(int count);

    private:
        // Konstruktor/Destruktor
        // Physics update functions
        void updateKerrPhotons(double a_spin);
        void updateTestPhotons(double h, double M);
        glm::vec3 getDopplerColor(const KerrState& s0, const KerrState& s1);
        // GPU compute update
        void initPhotonSSBOs();
        void uploadPhotonSSBOs();
        void dispatchPhotonCompute(float dt, float aSpin, float speedScale, float dopplerStrength);
        void readTrailsBackToCPU();

        // Rendering functions
        void drawKerrPhotons(const std::vector<Photon>& photons,
            GLuint shader, GLuint vao, GLuint vbo,
            const glm::mat4& proj, const glm::mat4& view,
            float zoom);
        void drawTestPhotons(GLuint shader, const glm::mat4 &projMat, const glm::mat4 &viewMat);
        
        // Photon initialization helpers
        void initializeKerrPhotons();
        void initializeTestPhotons();
        double calculatePhotonYPosition(int index, double yMin, double yMax) const;
        template <typename T>
        void trimPhotonTrail(std::vector<T>& trail) {
            if (trail.size() > trailLength) {
                const size_t removeCount = trail.size() - trailLength;
                trail.erase(trail.begin(), trail.begin() + removeCount);
            }
        }
        
        // Generic rendering helper
        template<typename PhotonType>
        void renderPhotonTrails(GLuint shader, const glm::mat4& projMat, const glm::mat4& viewMat,
                               const std::vector<PhotonType>& photons, const PhotonRenderConfig& config);
        void renderTrailSegments(int vertexCount, const PhotonRenderConfig& config);
        void outputDebugInfo();
    };
}
#endif // BLACKHOLESIM_ENGINE_HPP