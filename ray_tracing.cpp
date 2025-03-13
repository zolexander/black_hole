#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#include <cmath>
using namespace glm;

// vars
const int WIDTH = 800;
const int HEIGHT = 600;
glm::vec3 cameraPos   = glm::vec3(0.0f, 0.0f,  1.0f);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp    = glm::vec3(0.0f, 1.0f,  0.0f);
float lastX = 400.0, lastY = 300.0;
float cameraYaw = -90;
float cameraPitch = 0.0;
float deltaTime = 0.0;
float lastFrame = 0.0;

// functions
GLFWwindow* StartGLU();
GLuint CreateShaderProgram();
GLuint setupQuad();
void renderScene(GLFWwindow* &window, GLuint &quadVAO, GLuint &texture, GLuint &shaderProgram, std::vector<unsigned char> &pixels);
GLuint loadTexture();

void UpdateCam(GLuint shaderProgram, glm::vec3 cameraPos);
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);

// structures and classes :D
struct Ray{
    vec3 direction;
    vec3 origin;
    Ray(vec3 direction, vec3 origin) : direction(direction), origin(origin){}
};
struct Object{
    vec3 centre;
    float radius;
    vec3 color;
    Object(vec3 centre, float radius, vec3 color) : centre(centre), radius(radius), color(color){}
    bool intersect(const Ray &ray, float& t){
        vec3 oc = ray.origin - centre;
        double a = dot(ray.direction, ray.direction);
        double b = 2.0 * dot(oc, ray.direction);
        double c = dot(oc, oc) - radius * radius;

        double discriminant = b*b - 4 * a * c;
        if (discriminant < 0) return false;
        
        float temp = (-b - sqrt(discriminant)) / (2.0f*a);
        if (temp < 0) {
            temp = (-b + sqrt(discriminant)) / (2.0f*a);
            if (temp < 0) return false;
        }
        t = temp;
        return true;
    }
    
    vec3 getNormal(const glm::vec3& point) const {
        return normalize(point - centre);
    }
};
class Scene{
public:
    std::vector<Object> objs;
    vec3 lightPos;
    Scene() : lightPos(5.0f, 5.0f, 5.0f) {}

    vec3 trace(const Ray& ray) const{
        float closest = INFINITY;
        const Object* hitObject = nullptr;
        for(const auto& obj : objs){
            float t; // distance to intersection
            // Create non-const copy to call intersect
            Object mutableObj = obj;
            if(mutableObj.intersect(ray, t)) {
                if(t < closest) {
                    closest = t;
                    hitObject = &obj;
                }
            }
        };
        if (hitObject) {
            vec3 hitPoint = ray.origin + ray.direction * closest;
            vec3 normal = hitObject->getNormal(hitPoint);
            vec3 lightDir = normalize(lightPos - hitPoint);

            float diff = std::max(glm::dot(normal, lightDir), 0.0f);
            Ray shadowRay(lightDir, hitPoint + normal * 0.001f);
            bool inShadow = false;
            for (const auto& obj : objs) {
                float t;
                Object mutableObj = obj; // Create non-const copy
                if(mutableObj.intersect(shadowRay, t)) {
                    inShadow = true;
                    break;
                }
            }
            
            glm::vec3 color = hitObject->color; 
            float ambient = 0.1f; 
            
            if (inShadow) {
                return color * ambient;
            }
            
            return color * (ambient + diff * 0.9f);
        }
        return vec3(0.0f);
    }
};

// --------- main --------- //
int main() {
    // setup
    GLFWwindow* window = StartGLU();
    Scene scene;
    GLuint shaderProgram = CreateShaderProgram(); // compile shader program
    GLuint quadVAO = setupQuad(); // create quad background
    GLuint texture = loadTexture();
    glfwSetKeyCallback(window, keyCallback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    std::vector<unsigned char> pixels(WIDTH * HEIGHT * 3);
    scene.objs = {
        Object(glm::vec3(-1.0f, 0.0f, -5.0f), 1.0f, vec3(1.0f, 0.0, 0.0)),
        Object(glm::vec3(2.0f, 0.0f, -6.0f), 1.0f, vec3(0.0f, 1.0, 1.0))
    };
    
    while(!glfwWindowShouldClose(window)){
        glClear(GL_COLOR_BUFFER_BIT);
        UpdateCam(shaderProgram, cameraPos);

        // ray tracing
        for(int y = 0; y < HEIGHT; ++y){
            for(int x = 0; x < WIDTH; ++x){
                float u = float(x) / float(WIDTH); // normalize 0-1
                float v = float(y) / float(HEIGHT);  // normalize 0-1
                vec3 cameraRight = normalize(cross(cameraFront, cameraUp));

                float aspectRatio = float(WIDTH) / float(HEIGHT);
                float fov = 45.0f;
                float halfHeight = tan(radians(fov / 2.0f));
                float halfWidth = aspectRatio * halfHeight;

                vec3 direction = normalize( 
                    cameraFront
                    + (2.0f * u - 1.0f) * halfWidth * cameraRight
                    + (1.0f - 2.0f * v) * halfHeight * cameraUp
                );

                Ray ray(glm::normalize(direction), cameraPos);
                vec3 color = scene.trace(ray);
            
                int index = (y * WIDTH + x) * 3;
                pixels[index + 0] = static_cast<unsigned char>(color.r * 255);
                pixels[index + 1] = static_cast<unsigned char>(color.g * 255);
                pixels[index + 2] = static_cast<unsigned char>(color.b * 255);

            }   
        }

        // actualize scene
        renderScene(window, quadVAO, texture, shaderProgram, pixels);

    }
    glfwTerminate();
}

// function declarations
GLFWwindow* StartGLU() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return nullptr;
    }

    GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "RAY_TRACING", nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return nullptr;
    }
    glfwMakeContextCurrent(window);

    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW." << std::endl;
        glfwTerminate();
        return nullptr;
    }

    glViewport(0, 0, WIDTH, HEIGHT);
    return window;
}
GLuint CreateShaderProgram(){

    const char* vertexShaderSource = R"(
    #version 330 core
    layout (location = 0) in vec3 aPos;
    layout (location = 1) in vec2 aTexCoord;
    out vec2 TexCoord;
    void main() {
        gl_Position = vec4(aPos, 1.0);
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

    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, nullptr);
    glCompileShader(vertexShader);

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
}
GLuint setupQuad() {
    float quadVertices[] = {
        // positions   // texCoords
        -1.0f,  1.0f,  0.0f, 1.0f,
        -1.0f, -1.0f,  0.0f, 0.0f,
         1.0f, -1.0f,  1.0f, 0.0f,
        -1.0f,  1.0f,  0.0f, 1.0f,
         1.0f, -1.0f,  1.0f, 0.0f,
         1.0f,  1.0f,  1.0f, 1.0f
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

    return VAO;
}
GLuint loadTexture(){
    GLuint texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    return texture;
};
void renderScene(GLFWwindow* &window, GLuint &quadVAO, GLuint &texture, GLuint &shaderProgram, std::vector<unsigned char> &pixels){
        // Update texture with ray traced result
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, WIDTH, HEIGHT, 0, GL_RGB, 
                 GL_UNSIGNED_BYTE, pixels.data());
    glUseProgram(shaderProgram);
    glBindVertexArray(quadVAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glfwSwapBuffers(window);
    glfwPollEvents();
}

void UpdateCam(GLuint shaderProgram, glm::vec3 cameraPos) {
    // Calculate current frame time for smooth movement
    float currentFrame = glfwGetTime();
    deltaTime = currentFrame - lastFrame;
    lastFrame = currentFrame;
}
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    float cameraSpeed = 1.0f * deltaTime;
    
    if (key == GLFW_KEY_W && (action == GLFW_PRESS || action == GLFW_REPEAT))
        cameraPos += cameraSpeed * cameraFront;
    if (key == GLFW_KEY_S && (action == GLFW_PRESS || action == GLFW_REPEAT))
        cameraPos -= cameraSpeed * cameraFront;
    if (key == GLFW_KEY_A && (action == GLFW_PRESS || action == GLFW_REPEAT))
        cameraPos -= normalize(cross(cameraFront, cameraUp)) * cameraSpeed;
    if (key == GLFW_KEY_D && (action == GLFW_PRESS || action == GLFW_REPEAT))
        cameraPos += normalize(cross(cameraFront, cameraUp)) * cameraSpeed;
    if (key == GLFW_KEY_Q && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}
void mouse_callback(GLFWwindow* window, double xpos, double ypos) {

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; 
    lastX = xpos;
    lastY = ypos;

    float sensitivity = 0.1f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    cameraYaw += xoffset;
    cameraPitch += -yoffset;

    if(cameraPitch > 89.0f) cameraPitch = 89.0f;
    if(cameraPitch < -89.0f) cameraPitch = -89.0f;

    glm::vec3 front;
    front.x = cos(glm::radians(cameraYaw)) * cos(glm::radians(cameraPitch));
    front.y = sin(glm::radians(cameraPitch));
    front.z = sin(glm::radians(cameraYaw)) * cos(glm::radians(cameraPitch));
    cameraFront = glm::normalize(front);
}
