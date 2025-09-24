#include "engine.hpp"
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <chrono>
#include "engine.hpp"
#include <kerrstate.hpp>
#include <kerrintegrate.hpp>
#include <kerr_inline.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
int main() {
    BlackholeSim::Engine engine;
    if (!engine.initGL("Photon Trails Demo", 1200, 800)) {
        return -1;
    }
    engine.run();
    return 0;
}