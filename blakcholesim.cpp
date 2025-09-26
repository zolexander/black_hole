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
#include <cmath>
// glm
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <blackholesim.hpp>

BlackholeSim::Photon::Photon(double L_ = 0.0) : L(L_), h(0.01), alive(true)
{
    s = {4.0, M_PI};
    trail.reserve(256);
}
void BlackholeSim::Photon::reset(double L_, const KerrState &startState)
{
    L = L_;
    s = startState;
    h = 0.01;
    alive = true;
    trail.clear();
    trail.reserve(256);
}
BlackholeSim::TestPhoton::TestPhoton(double x0 = -20.0, double y0 = 0.0)
{
    s = {x0, y0, 1.0, 0.0};
    alive = true;
    trail.reserve(256);
}
void BlackholeSim::TestPhoton::reset(double x0, double y0)
{
    s.x = x0;
    s.y = y0;
    s.vx = 1.0;
    s.vy = 0.0;
    alive = true;
    trail.clear();
    trail.reserve(256);
}