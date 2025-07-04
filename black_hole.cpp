#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <chrono>
#include <fstream>
#include <sstream>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace glm;
using namespace std;
using Clock = std::chrono::high_resolution_clock;

// VARS
double lastPrintTime = 0.0;
int    framesCount   = 0;
double c = 299792458.0;
double G = 6.67430e-11;
bool useGeodesics = true;

struct Camera {
    vec3 target = vec3(0.0f);           // Point the camera looks at
    float radius = 6.34194e10f;         // Distance from target
    float minRadius = 1e10f, maxRadius = 1e12f;

    float azimuth = 0.0f;               // Horizontal rotation
    float elevation = M_PI / 2.0f;      // Vertical rotation

    float orbitSpeed = 0.01f;
    float panSpeed = 0.01f; // reduced to avoid excessive movement
    double zoomSpeed = 25e9f;

    bool dragging = false;
    bool panning = false;
    double lastX = 0.0, lastY = 0.0;

    // Calculate camera position in world space
    vec3 position() const {
        float clampedElevation = clamp(elevation, 0.01f, float(M_PI) - 0.01f);
        return target + vec3(
            radius * sin(clampedElevation) * cos(azimuth),
            radius * cos(clampedElevation),
            radius * sin(clampedElevation) * sin(azimuth)
        );
    }
    void update() {
        // You could add smoothing/caching here if needed
    }

    void processMouseMove(double x, double y) {
        float dx = float(x - lastX);
        float dy = float(y - lastY);

        if (dragging && panning) {
            // Pan: Shift + Left or Middle Mouse
            vec3 camPos = position();
            vec3 forward = normalize(target - camPos);
            vec3 worldUp = vec3(0, 1, 0);
            vec3 right = normalize(cross(forward, worldUp));
            vec3 up = normalize(cross(right, forward));

            // Inverted direction to match drag feel
            target -= right * dx * panSpeed * radius;
            target += up * dy * panSpeed * radius;
        } 
        else if (dragging && !panning) {
            // Orbit: Left mouse only
            azimuth   += dx * orbitSpeed;
            elevation -= dy * orbitSpeed;
            elevation = clamp(elevation, 0.01f, float(M_PI) - 0.01f);
        }

        lastX = x;
        lastY = y;
        update();
    }
    void processMouseButton(int button, int action, int mods, GLFWwindow* win) {
        if (button == GLFW_MOUSE_BUTTON_LEFT || button == GLFW_MOUSE_BUTTON_MIDDLE) {
            if (action == GLFW_PRESS) {
                dragging = true;
                panning = (button == GLFW_MOUSE_BUTTON_MIDDLE) || 
                         (button == GLFW_MOUSE_BUTTON_LEFT && (mods & GLFW_MOD_SHIFT));
                glfwGetCursorPos(win, &lastX, &lastY);
            } else if (action == GLFW_RELEASE) {
                dragging = false;
                panning = false;
            }
        }
    }
    void processScroll(double xoffset, double yoffset) {
        radius -= yoffset * zoomSpeed;
        radius = clamp(radius, minRadius, maxRadius);
        update();
    }
};
Camera camera;

struct Ray;
void rk4Step(Ray& ray, double dλ, double rs);
struct BlackHole {
    vec3 position;
    double mass;
    double radius;
    double r_s;

    BlackHole(vec3 pos, float m) : position(pos), mass(m) {r_s = 2.0 * G * mass / (c*c);}
    bool Intercept(float px, float py, float pz) const {
        double dx = double(px) - double(position.x);
        double dy = double(py) - double(position.y);
        double dz = double(pz) - double(position.z);
        double dist2 = dx * dx + dy * dy + dz * dz;
        return dist2 < r_s * r_s;
    }
};
BlackHole SagA(vec3(0.0f, 0.0f, 0.0f), 8.54e36); // Sagittarius A black hole
struct Ray{
    // -- cartesian coords -- //
    double x;   double y; double z;
    // -- polar coords -- //
    double r;   double phi; double theta;
    double dr;  double dphi; double dtheta;
    double E, L;             // conserved quantities

    Ray(vec3 pos, vec3 dir) : x(pos.x), y(pos.y), z(pos.z) {
        // Step 1: get spherical coords (r, theta, phi)
        r = sqrt(x*x + y*y + z*z);
        theta = acos(z / r);
        phi = atan2(y, x);

        // Step 2: seed velocities (dr, dtheta, dphi)
        // Convert direction to spherical basis
        double dx = dir.x, dy = dir.y, dz = dir.z;
        dr     = sin(theta)*cos(phi)*dx + sin(theta)*sin(phi)*dy + cos(theta)*dz;
        dtheta = cos(theta)*cos(phi)*dx + cos(theta)*sin(phi)*dy - sin(theta)*dz;
        dtheta /= r;
        dphi   = -sin(phi)*dx + cos(phi)*dy;
        dphi  /= (r * sin(theta));

        // Step 3: store conserved quantities
        L = r * r * sin(theta) * dphi;
        double f = 1.0 - SagA.r_s / r;
        double dt_dλ = sqrt((dr*dr)/f + r*r*dtheta*dtheta + r*r*sin(theta)*sin(theta)*dphi*dphi);
        E = f * dt_dλ;
    }
    void step(double dλ, double rs) {
        if (r <= rs) return;
        rk4Step(*this, dλ, rs);
        // convert back to cartesian
        this->x = r * sin(theta) * cos(phi);
        this->y = r * sin(theta) * sin(phi);
        this->z = r * cos(theta);
    }
};

struct Engine {
    // -- Quad & Texture render -- //
    GLFWwindow* window;
    GLuint quadVAO;
    GLuint texture;
    GLuint shaderProgram;
    GLuint computeProgram = 0;
    GLuint cameraUBO = 0;
    GLuint diskUBO = 0;
    int WIDTH = 800;  // Window width
    int HEIGHT = 600; // Window height
    int COMPUTE_WIDTH = 400;   // Compute resolution width
    int COMPUTE_HEIGHT = 300;  // Compute resolution height
    float width = 100000000000.0f; // Width of the viewport in meters
    float height = 75000000000.0f; // Height of the viewport in meters
    
    Engine() {
        if (!glfwInit()) {
            cerr << "GLFW init failed\n";
            exit(EXIT_FAILURE);
        }
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        window = glfwCreateWindow(WIDTH, HEIGHT, "Black Hole", nullptr, nullptr);
        if (!window) {
            cerr << "Failed to create GLFW window\n";
            glfwTerminate();
            exit(EXIT_FAILURE);
        }
        glfwMakeContextCurrent(window);
        glewExperimental = GL_TRUE;
        GLenum glewErr = glewInit();
        if (glewErr != GLEW_OK) {
            cerr << "Failed to initialize GLEW: "
                << (const char*)glewGetErrorString(glewErr)
                << "\n";
            glfwTerminate();
            exit(EXIT_FAILURE);
        }
        cout << "OpenGL " << glGetString(GL_VERSION) << "\n";
        this->shaderProgram = CreateShaderProgram();

        computeProgram = CreateComputeProgram("geodesic.comp");
        glGenBuffers(1, &cameraUBO);
        glBindBuffer(GL_UNIFORM_BUFFER, cameraUBO);
        glBufferData(GL_UNIFORM_BUFFER, 128, nullptr, GL_DYNAMIC_DRAW); // alloc ~128 bytes
        glBindBufferBase(GL_UNIFORM_BUFFER, 1, cameraUBO); // binding = 1 matches shader

        glGenBuffers(1, &diskUBO);
        glBindBuffer(GL_UNIFORM_BUFFER, diskUBO);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(float) * 4, nullptr, GL_DYNAMIC_DRAW); // 3 values + 1 padding
        glBindBufferBase(GL_UNIFORM_BUFFER, 2, diskUBO); // binding = 2 matches compute shader

        auto result = QuadVAO();
        this->quadVAO = result[0];
        this->texture = result[1];
    }
    GLuint CreateShaderProgram(){
        const char* vertexShaderSource = R"(
        #version 330 core
        layout (location = 0) in vec2 aPos;  // Changed to vec2
        layout (location = 1) in vec2 aTexCoord;
        out vec2 TexCoord;
        void main() {
            gl_Position = vec4(aPos, 0.0, 1.0);  // Explicit z=0
            TexCoord = aTexCoord;
        })";

        const char* fragmentShaderSource = R"(
        #version 330 core
        in vec2 TexCoord;
        out vec4 FragColor;
        uniform sampler2D screenTexture;
        void main() {
            FragColor = texture(screenTexture, TexCoord);
        })";

        // vertex shader
        GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertexShader, 1, &vertexShaderSource, nullptr);
        glCompileShader(vertexShader);

        // fragment shader
        GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragmentShader, 1, &fragmentShaderSource, nullptr);
        glCompileShader(fragmentShader);

        GLuint shaderProgram = glCreateProgram();
        glAttachShader(shaderProgram, vertexShader);
        glAttachShader(shaderProgram, fragmentShader);
        glLinkProgram(shaderProgram);

        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);

        return shaderProgram;
    };
    GLuint CreateComputeProgram(const char* path) {
        // 1) read GLSL source
        std::ifstream in(path);
        if(!in.is_open()) {
            std::cerr << "Failed to open compute shader: " << path << "\n";
            exit(EXIT_FAILURE);
        }
        std::stringstream ss;
        ss << in.rdbuf();
        std::string srcStr = ss.str();
        const char* src = srcStr.c_str();

        // 2) compile
        GLuint cs = glCreateShader(GL_COMPUTE_SHADER);
        glShaderSource(cs, 1, &src, nullptr);
        glCompileShader(cs);
        GLint ok; 
        glGetShaderiv(cs, GL_COMPILE_STATUS, &ok);
        if(!ok) {
            GLint logLen;
            glGetShaderiv(cs, GL_INFO_LOG_LENGTH, &logLen);
            std::vector<char> log(logLen);
            glGetShaderInfoLog(cs, logLen, nullptr, log.data());
            std::cerr << "Compute shader compile error:\n" << log.data() << "\n";
            exit(EXIT_FAILURE);
        }

        // 3) link
        GLuint prog = glCreateProgram();
        glAttachShader(prog, cs);
        glLinkProgram(prog);
        glGetProgramiv(prog, GL_LINK_STATUS, &ok);
        if(!ok) {
            GLint logLen;
            glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &logLen);
            std::vector<char> log(logLen);
            glGetProgramInfoLog(prog, logLen, nullptr, log.data());
            std::cerr << "Compute shader link error:\n" << log.data() << "\n";
            exit(EXIT_FAILURE);
        }

        glDeleteShader(cs);
        return prog;
    }
    void dispatchCompute(const Camera& cam) {
        // 1) bind your compute pipeline
        glUseProgram(computeProgram);

        uploadCameraUBO(cam);

        // disk
        float r1 = SagA.r_s * 2.2f;    // inner radius just outside the event horizon
        float r2 = SagA.r_s * 4.2f;   // outer radius of the disk
        float q  = 2.0f;               // emissivity falloff (brightness ∝ 1/r^q)
        float padding = 0.0f;          // padding for std140 alignment
        float diskData[4] = { r1, r2, q, padding };

        glBindBuffer(GL_UNIFORM_BUFFER, diskUBO);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(diskData), diskData);

        // 3) bind your render‐texture as image unit 0
        glBindImageTexture(0, texture, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA8);

        // 4) launch the compute grid (16×16 workgroups)
        GLuint groupsX = (GLuint)std::ceil(static_cast<float>(COMPUTE_WIDTH) / 16.0f);
        GLuint groupsY = (GLuint)std::ceil(static_cast<float>(COMPUTE_HEIGHT) / 16.0f);
        glDispatchCompute(groupsX, groupsY, 1);

        // 5) make sure writes are visible to the rendering pipeline
        glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
    }
    void uploadCameraUBO(const Camera& cam) {
        struct UBOData {
            vec3 pos; float _pad0;
            vec3 right; float _pad1;
            vec3 up; float _pad2;
            vec3 forward; float _pad3;
            float tanHalfFov;
            float aspect;
            int useGeodesics;
            int _pad4;
        } data;
        vec3 fwd = normalize(cam.target - cam.position());
        vec3 right = normalize(cross(fwd, vec3(0,1,0)));
        vec3 up = cross(right, fwd);

        data.pos = cam.position();     data.right = right;
        data.up = up;           data.forward = fwd;
        data.tanHalfFov = tan(radians(60.0f * 0.5f));
        data.aspect = float(WIDTH) / float(HEIGHT); // Window aspect ratio
        data.useGeodesics = useGeodesics;

        glBindBuffer(GL_UNIFORM_BUFFER, cameraUBO);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(UBOData), &data);
    }
    vector<GLuint> QuadVAO(){
        float quadVertices[] = {
            // positions   // texCoords
            -1.0f,  1.0f,  0.0f, 1.0f,  // top left
            -1.0f, -1.0f,  0.0f, 0.0f,  // bottom left
            1.0f, -1.0f,  1.0f, 0.0f,  // bottom right

            -1.0f,  1.0f,  0.0f, 1.0f,  // top left
            1.0f, -1.0f,  1.0f, 0.0f,  // bottom right
            1.0f,  1.0f,  1.0f, 1.0f   // top right
        };
        
        GLuint VAO, VBO;
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);

        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);

        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
        glEnableVertexAttribArray(1);

        GLuint texture;
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, COMPUTE_WIDTH, COMPUTE_HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
        vector<GLuint> VAOtexture = {VAO, texture};
        return VAOtexture;
    }
    void renderScene() {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glUseProgram(shaderProgram);
        glBindVertexArray(quadVAO);
        // make sure your fragment shader samples from texture unit 0:
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texture);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        glfwSwapBuffers(window);
        glfwPollEvents();
    };
    static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
        if (action == GLFW_PRESS) {
            if (key == GLFW_KEY_G) {
                useGeodesics = !useGeodesics;
                cout << "Geodesics: " << (useGeodesics ? "ON\n" : "OFF\n");
            }
        }
    }
};
Engine engine;
void setupCameraCallbacks(GLFWwindow* window) {
    glfwSetWindowUserPointer(window, &camera); // So callbacks can access the camera

    glfwSetMouseButtonCallback(window, [](GLFWwindow* win, int button, int action, int mods) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processMouseButton(button, action, mods, win);
    });

    glfwSetCursorPosCallback(window, [](GLFWwindow* win, double x, double y) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processMouseMove(x, y);
    });

    glfwSetScrollCallback(window, [](GLFWwindow* win, double xoffset, double yoffset) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processScroll(xoffset, yoffset);
    });
    glfwSetKeyCallback(window, [](GLFWwindow* win, int key, int sc, int action, int mods){
        Engine::keyCallback(win, key, sc, action, mods);
    });
}

// -- MAIN -- //
int main() {
    setupCameraCallbacks(engine.window);
    vector<unsigned char> pixels(engine.WIDTH * engine.HEIGHT * 3);

    auto t0 = Clock::now();
    lastPrintTime = chrono::duration<double>(t0.time_since_epoch()).count();

    while (!glfwWindowShouldClose(engine.window)) {
        engine.dispatchCompute(camera);
        engine.renderScene();

        // 2) FPS counting
        framesCount++;
        auto t1 = Clock::now();
        double now = chrono::duration<double>(t1.time_since_epoch()).count();
        if (now - lastPrintTime >= 1.0) {
            cout << "FPS: " << framesCount / (now - lastPrintTime) << "\n";
            cout << "Camera Position: " << camera.position().x << ", "
                 << camera.position().y << ", " << camera.position().z<<endl;
            cout << "radius: " << camera.radius << "\n";
            framesCount   = 0;
            lastPrintTime = now;
        }
    }

    glfwDestroyWindow(engine.window);
    glfwTerminate();
    return 0;
}
