#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
using namespace glm;

GLFWwindow* StartGLU();
struct Ray{
    vec3 direction;
    vec3 origin;
    Ray(vec3 direction, vec3 origin) : direction(direction), origin(origin) {}
};
struct Object{
    float radius;
    vec3 centre;
    Object(float radius, vec3 centre) : centre(centre), radius(radius){}
};

int main() {
    GLFWwindow* window = StartGLU();

    while(!glfwWindowShouldClose(window)){

        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    glfwTerminate();
}

GLFWwindow* StartGLU() {
    if (!glfwInit()) {
        std::cout << "Failed to initialize GLFW, panic" << std::endl;
        return nullptr;
    }
    GLFWwindow* window = glfwCreateWindow(800, 600, "RAY_TRACING", NULL, NULL);
    if (!window) {
        std::cerr << "Failed to create GLFW window, PANIC PANIC PANIC!" << std::endl;
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

    glEnable(GL_DEPTH_TEST);
    glViewport(0, 0, 800, 600);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); // Standard blending for transparency

    return window;
}





