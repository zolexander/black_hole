#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace glm;
using namespace std;

// --- VARS --- //
float G = 6.67430 * pow(10, -11);
float solarMass = 1.989 * pow(10, 30); // kg
float M = solarMass * 4297000.0f; // kg
float c = 299792458.0f; // m/s

vec2 blackHolePos = vec2(400.0f, 300.0f);

// --- CONVERSIONS --- //
float MtoPx = c / 1000.0f; // pixels per km - lights travels 1kpx/s
float PxToM = 1.0f / MtoPx; // pixels to meters

struct Ray;
void rk4Step(Ray& ray, double dλ, double rs);

// --- STRUCTS --- //
struct Ray {
    // -- cartesian coords -- //
    double x;   double y;
    // -- polar coords -- //
    double r;   double phi;
    double dr;  double dphi;
    std::vector<vec2> trail;  // ← store past positions

    Ray(vec2 pos, vec2 vel) : x(pos[0]) {
        // step 1) store pos & dist : 
        x = pos.x; y = pos.y;
        double dx_m = (x - blackHolePos.x) * MtoPx;
        double dy_m = (y - blackHolePos.y) * MtoPx;

        // step 2) get polar coords (r, phi) :
        r = sqrt(dx_m * dx_m + dy_m * dy_m); // m
        phi = atan2(dy_m, dx_m); // radians

        // step 3) seed velocities :
        dr = vel.x * cos(phi) + vel.y * sin(phi); // m/s
        dphi  = ( -vel.x * sin(phi) + vel.y * cos(phi) ) / r;

        // step 4) start trail : 
        trail.push_back({x, y});
    }
    
    void step(double dλ, double rs) {
        rk4Step(*this, dλ, rs);  // evolve physics (r,φ,dr,dφ)

        // Convert back to screen space
        x = blackHolePos.x + r * cos(phi) * PxToM;
        y = blackHolePos.y + r * sin(phi) * PxToM;

        trail.push_back({x, y});
    }
};
struct physics {
    float r_s = 2 * G * M / (c * c); // m
    float pxWidth = 15.0f;           // width of black hole in pixels
    
     // Schwarzschild radius
    float meterToPx = pxWidth / r_s;    // 1 meter to pixels
    float pxToMeter = r_s / pxWidth;    // pixels to meters

    float rs_m() {
        return r_s * meterToPx; // convert to pixels
    }
};

struct Engine {
    GLFWwindow* window;
    int WIDTH = 800;
    int HEIGHT = 600;

    std::vector<Ray>* raysPtr = nullptr;

    Engine() {
        window = StartGLFW();
        glfwSetWindowUserPointer(this->window, this);
        glfwSetMouseButtonCallback(this->window, this->mouseButtonCallback);
    }

    GLFWwindow* StartGLFW() {
        if (!glfwInit()) {
            std::cerr << "glfw failed init, PANIC PANIC!" << std::endl;
            return nullptr;
        }
    
        GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "ray tracer", NULL, NULL);
        glfwMakeContextCurrent(window);
    
        glewExperimental = GL_TRUE;
        if (glewInit() != GLEW_OK) {
            std::cerr << "Failed to initialize GLEW." << std::endl;
            glfwTerminate();
            return nullptr;
        }
    
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glViewport(0, 0, WIDTH, HEIGHT);
    
        // ←––– ADD THIS: set up a 0…WIDTH, 0…HEIGHT orthographic projection
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(0.0, WIDTH,    // left, right
                0.0, HEIGHT,   // bottom, top
               -1.0, 1.0);     // near, far
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    
        return window;
    }
    
    void renderScreen(float r, std::vector<Ray>& rays) {
        this->draw(rays);
        this->drawCircle(blackHolePos.x, blackHolePos.y, r, 100);
        glfwSwapBuffers(window);
        glfwPollEvents();              
    }
    void drawCircle(float cx, float cy, float r, int segments) {
        glColor3f(1.0f, 0.0f, 0.0f);
        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(cx, cy); // center
        for (int i = 0; i <= segments; ++i) {
            float angle = 2.0f * M_PI * i / segments;
            float x = r * cos(angle);
            float y = r * sin(angle);
            glVertex2f(cx + x, cy + y);
        }
        glEnd();
    }
    void draw(const std::vector<Ray>& rays) {
        // draw current ray positions as points
        glPointSize(2.0f);
        glColor3f(0.0f, 0.0f, 0.0f);
        glBegin(GL_POINTS);
          for (const auto& ray : rays) {
              glVertex2f(ray.x, ray.y);
          }
        glEnd();
    
        // turn on blending for the trails
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glLineWidth(1.0f);
    
        // draw each trail with fading alpha
        for (const auto& ray : rays) {
            size_t N = ray.trail.size();
            if (N < 2) continue;
    
            glBegin(GL_LINE_STRIP);
            for (size_t i = 0; i < N; ++i) {
                // older points (i=0) get alpha≈0, newer get alpha≈1
                float alpha = float(i) / float(N - 1);
                glColor4f(1.0f, 1.0f, 1.0f, std::max(alpha, 1.05f));
                glVertex2f(ray.trail[i].x, ray.trail[i].y);
            }
            glEnd();
        }
    
        glDisable(GL_BLEND);
    }
    static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
        if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
            Engine* engine = static_cast<Engine*>(glfwGetWindowUserPointer(window));
            if (engine->raysPtr) {
                for (int i = 0; i < 100; ++i) {
                    float y = 100 + i * 4;
                    engine->raysPtr->push_back({ vec2(800.0f, y), vec2(-c / 100.0f, 0.0f) });
                }
            }
        }
    }
};

void geodesicRHS(const Ray& ray, double rhs[4], double rs) {
    double r    = ray.r;
    double dr   = ray.dr;
    double dphi = ray.dphi;

    double f = 1.0 - rs / r;

    rhs[0] = dr;
    rhs[1] = dphi;
    rhs[2] = -(rs / (2.0 * r * r)) * ( (dr * dr) / f + f * r * r * dphi * dphi );
    rhs[3] = -2.0 * dr * dphi / r;
}
void addState(const double a[4], const double b[4], double factor, double out[4]) {
    for (int i = 0; i < 4; i++)
        out[i] = a[i] + b[i] * factor;
}
void rk4Step(Ray& ray, double dλ, double rs) {
    double y0[4] = { ray.r, ray.phi, ray.dr, ray.dphi };
    double k1[4], k2[4], k3[4], k4[4], temp[4];

    // k1
    geodesicRHS(ray, k1, rs);

    // k2
    addState(y0, k1, dλ/2.0, temp);
    Ray r2 = ray;
    r2.r = temp[0]; r2.phi = temp[1]; r2.dr = temp[2]; r2.dphi = temp[3];
    geodesicRHS(r2, k2, rs);

    // k3
    addState(y0, k2, dλ/2.0, temp);
    Ray r3 = ray;
    r3.r = temp[0]; r3.phi = temp[1]; r3.dr = temp[2]; r3.dphi = temp[3];
    geodesicRHS(r3, k3, rs);

    // k4
    addState(y0, k3, dλ, temp);
    Ray r4 = ray;
    r4.r = temp[0]; r4.phi = temp[1]; r4.dr = temp[2]; r4.dphi = temp[3];
    geodesicRHS(r4, k4, rs);

    // Final update
    ray.r   += (dλ/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
    ray.phi   += (dλ/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    ray.dr  += (dλ/6.0) * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);
    ray.dphi  += (dλ/6.0) * (k1[3] + 2*k2[3] + 2*k3[3] + k4[3]);
}


// --- MAIN --- //
int main() {
    Engine engine;
    physics phys;
    std::vector<Ray> rays;

    // Init rays
    for (int i = 0; i < 150; ++i) {
        float y = 50 + i * 4;
        rays.push_back({ vec2(800.0f, y), vec2(-c / 100.0f, 0.0f) });
    }
    //rays.push_back({ vec2(800.0f, 300.0f), vec2(-c, 0.0f) });

    // point engine to rays vector
    engine.raysPtr = &rays;
    double rs_m = phys.r_s;

    // --- LOOP --- //
    while (!glfwWindowShouldClose(engine.window)) {
        glClear(GL_COLOR_BUFFER_BIT);
        glLoadIdentity();

        // Update rays
        float dt = 0.1f;
        for (auto& ray : rays) {

            // - Check Event Horizon - //
            if (length(vec2(ray.x, ray.y) - blackHolePos) < phys.rs_m()) continue;

            ray.step(dt, phys.rs_m()); // evolve physics (r,φ,dr,dφ)

        }


        // Render scene
        engine.renderScreen(phys.rs_m(), rays);
    }

    return 0;
}
