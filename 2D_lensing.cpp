#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"
#define _USE_MATH_DEFINES
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace glm;
using namespace std;

// --- Physikalische Konstanten --- //
const double c = 299792458.0;
const double G = 6.67430e-11;

// --- Vorwärtsdeklarationen --- //
struct Ray;
void rk4Step(Ray& ray, double dLambda, double schwarzschildRadius);
void rk4StepKerr(Ray& ray, double dLambda, double rs, double a); // <--- Diese Zeile ergänzen

// --- Engine für Fenster und Navigation --- //
struct Engine {
    GLFWwindow* window;
    int pixelWidth = 800;
    int pixelHeight = 600;
    float physWidth = 1e11f;   // Meter
    float physHeight = 7.5e10f; // Meter

    // Navigation
    float offsetX = 0.0f, offsetY = 0.0f;
    float zoom = 1.0f;
    bool middleMousePressed = false;
    double lastMouseX = 0.0, lastMouseY = 0;

    Engine() {
        if (!glfwInit()) {
            cerr << "GLFW konnte nicht initialisiert werden." << endl;
            exit(EXIT_FAILURE);
        }
        window = glfwCreateWindow(pixelWidth, pixelHeight, "Black Hole Simulation", NULL, NULL);
        if (!window) {
            cerr << "GLFW Fenster konnte nicht erstellt werden." << endl;
            glfwTerminate();
            exit(EXIT_FAILURE);
        }
        glfwMakeContextCurrent(window);
        if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
            cerr << "GLAD konnte nicht initialisiert werden." << endl;
            glfwDestroyWindow(window);
            glfwTerminate();
            exit(EXIT_FAILURE);
        }
        glViewport(0, 0, pixelWidth, pixelHeight);
    }
    void setupImGui() {
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImGuiIO& io = ImGui::GetIO(); (void)io;
        ImGui::StyleColorsDark();
        ImGui_ImplGlfw_InitForOpenGL(window, true);
        ImGui_ImplOpenGL3_Init("#version 330");
    }
    void setupViewport() {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        double left   = -physWidth + offsetX;
        double right  =  physWidth + offsetX;
        double bottom = -physHeight + offsetY;
        double top    =  physHeight + offsetY;
        glOrtho(left, right, bottom, top, -1.0, 1.0);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    }
};
Engine engine;

// --- Schwarzes Loch --- //
struct BlackHole {
    vec3 position;
    double mass;
    double r_s; // Schwarzschild-Radius

    BlackHole(vec3 pos, double m) : position(pos), mass(m) {
        r_s = 2.0 * G * mass / (c * c);
    }

    void draw() {
        glBegin(GL_TRIANGLE_FAN);
        glColor3f(1.0f, 0.0f, 0.0f); // Rot
        glVertex2f(0.0f, 0.0f);      // Zentrum
        for (int i = 0; i <= 100; i++) {
            float angle = 2.0f * M_PI * i / 100;
            float x = r_s * cos(angle);
            float y = r_s * sin(angle);
            glVertex2f(x, y);
        }
        glEnd();
    }
};
BlackHole SagA(vec3(0.0f, 0.0f, 0.0f), 8.54e36);

// --- Lichtstrahl (Ray) --- //
struct Ray {
    // Kartesische Koordinaten
    double x, y;
    // Polarkoordinaten
    double r, phi;
    double dr, dphi;
    vector<vec2> trail; // Spur der Punkte
    double E, L;        // Erhaltungsgrößen

    Ray(vec2 pos, vec2 dir) : x(pos.x), y(pos.y) {
        // 1. Polarkoordinaten berechnen
        r = sqrt(x * x + y * y);
        phi = atan2(y, x);

        // 2. Anfangsgeschwindigkeiten in Polarkoordinaten
        dr = dir.x * cos(phi) + dir.y * sin(phi);
        dphi = (-dir.x * sin(phi) + dir.y * cos(phi)) / r;

        // 3. Erhaltungsgrößen berechnen
        L = r * r * dphi;
        double f = 1.0 - SagA.r_s / r;
        double dt_dLambda = sqrt((dr * dr) / (f * f) + (r * r * dphi * dphi) / f);
        E = f * dt_dLambda;

        // 4. Startpunkt zur Spur hinzufügen
        trail.push_back({x, y});
    }

    void drawTrailAndPoint(const std::vector<Ray>& rays) {
        // Aktuelle Positionen als Punkte zeichnen
        glPointSize(2.0f);
        glColor3f(1.0f, 0.0f, 0.0f);
        glBegin(GL_POINTS);
        for (const auto& ray : rays) {
            glVertex2f(ray.x, ray.y);
        }
        glEnd();

        // Trails mit Transparenz zeichnen
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glLineWidth(1.0f);

        for (const auto& ray : rays) {
            size_t N = ray.trail.size();
            if (N < 2) continue;
            glBegin(GL_LINE_STRIP);
            for (size_t i = 0; i < N; ++i) {
                float alpha = float(i) / float(N - 1);
                glColor4f(1.0f, 1.0f, 1.0f, std::max(alpha, 0.05f));
                glVertex2f(ray.trail[i].x, ray.trail[i].y);
            }
            glEnd();
        }
        glDisable(GL_BLEND);
    }

    // ...im Ray struct...
    void step(double dLambda, double schwarzschildRadius, double kerr_a = 0.0, bool useKerr = false) {
        if (r <= schwarzschildRadius) return;
        if (useKerr && kerr_a > 0.0) {
            rk4StepKerr(*this, dLambda, schwarzschildRadius, kerr_a);
        } else {
        rk4Step(*this, dLambda, schwarzschildRadius);
        }
        x = r * cos(phi);
        y = r * sin(phi);
        trail.push_back({float(x), float(y)});
    }
    };
vector<Ray> rays;
// Kerr-Äquator-Gleichung (vereinfachte Näherung für kleine a)
void geodesicRHS_Kerr(const Ray& ray, double rhs[4], double rs, double a) {
    // rs = 2GM/c^2, a = Kerr-Parameter (J/Mc)
    double r = ray.r;
    double dr = ray.dr;
    double dphi = ray.dphi;
    double u = 1.0 / r;
    double phi = ray.phi;

    // Näherung: d^2u/dphi^2 + u = 3Mu^2 - 2aMu^3
    // Wir schreiben das als System erster Ordnung:
    // y[0] = u, y[1] = du/dphi
    // dy[0]/dphi = y[1]
    // dy[1]/dphi = -y[0] + 3*M*y[0]^2 - 2*a*M*y[0]^3

    double M = rs / 2.0; // rs = 2GM/c^2 -> M = GM/c^2

    rhs[0] = dr; // dr/dlambda
    rhs[1] = dphi; // dphi/dlambda

    // du/dphi = dr/dlambda / dphi/dlambda
    double du_dphi = (dr / dphi) / (r * r);

    // d^2u/dphi^2
    double d2u_dphi2 = -u + 3.0 * M * u * u - 2.0 * a * M * u * u * u;

    // dr/dlambda = du/dphi * dphi/dlambda * r^2
    rhs[2] = d2u_dphi2 * dphi * r * r;
    // dphi bleibt gleich
    rhs[3] = 0.0;
}
// --- Geodätengleichungen (Schwarzschild) --- //
void geodesicRHS(const Ray& ray, double rhs[4], double rs) {
    double r = ray.r;
    double dr = ray.dr;
    double dphi = ray.dphi;
    double E = ray.E;
    double f = 1.0 - rs / r;

    rhs[0] = dr;
    rhs[1] = dphi;

    double dt_dLambda = E / f;
    rhs[2] = - (rs / (2 * r * r)) * f * (dt_dLambda * dt_dLambda)
             + (rs / (2 * r * r * f)) * (dr * dr)
             + (r - rs) * (dphi * dphi);

    rhs[3] = -2.0 * dr * dphi / r;
}

void addState(const double a[4], const double b[4], double factor, double out[4]) {
    for (int i = 0; i < 4; i++)
        out[i] = a[i] + b[i] * factor;
}

void rk4Step(Ray& ray, double dLambda, double rs) {
    double y0[4] = { ray.r, ray.phi, ray.dr, ray.dphi };
    double k1[4], k2[4], k3[4], k4[4], temp[4];

    geodesicRHS(ray, k1, rs);
    addState(y0, k1, dLambda / 2.0, temp);
    Ray r2 = ray; r2.r = temp[0]; r2.phi = temp[1]; r2.dr = temp[2]; r2.dphi = temp[3];
    geodesicRHS(r2, k2, rs);

    addState(y0, k2, dLambda / 2.0, temp);
    Ray r3 = ray; r3.r = temp[0]; r3.phi = temp[1]; r3.dr = temp[2]; r3.dphi = temp[3];
    geodesicRHS(r3, k3, rs);

    addState(y0, k3, dLambda, temp);
    Ray r4 = ray; r4.r = temp[0]; r4.phi = temp[1]; r4.dr = temp[2]; r4.dphi = temp[3];
    geodesicRHS(r4, k4, rs);

    ray.r    += (dLambda / 6.0) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    ray.phi  += (dLambda / 6.0) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    ray.dr   += (dLambda / 6.0) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
    ray.dphi += (dLambda / 6.0) * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]);
}
void rk4StepKerr(Ray& ray, double dLambda, double rs, double a) {
    double y0[4] = { ray.r, ray.phi, ray.dr, ray.dphi };
    double k1[4], k2[4], k3[4], k4[4], temp[4];

    geodesicRHS_Kerr(ray, k1, rs, a);
    addState(y0, k1, dLambda / 2.0, temp);
    Ray r2 = ray; r2.r = temp[0]; r2.phi = temp[1]; r2.dr = temp[2]; r2.dphi = temp[3];
    geodesicRHS_Kerr(r2, k2, rs, a);

    addState(y0, k2, dLambda / 2.0, temp);
    Ray r3 = ray; r3.r = temp[0]; r3.phi = temp[1]; r3.dr = temp[2]; r3.dphi = temp[3];
    geodesicRHS_Kerr(r3, k3, rs, a);

    addState(y0, k3, dLambda, temp);
    Ray r4 = ray; r4.r = temp[0]; r4.phi = temp[1]; r4.dr = temp[2]; r4.dphi = temp[3];
    geodesicRHS_Kerr(r4, k4, rs, a);

    ray.r    += (dLambda / 6.0) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    ray.phi  += (dLambda / 6.0) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    ray.dr   += (dLambda / 6.0) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
    ray.dphi += (dLambda / 6.0) * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]);
}
// --- Hauptprogramm --- //
int main() {
    // Strahlen initialisieren
    for (float y = -engine.physHeight; y < engine.physHeight; y += 1e10) {
        rays.push_back(Ray(vec2(-engine.physWidth, y), vec2(c, 0.0)));
    }
    engine.setupImGui();
    
    // Hauptloop
    while (!glfwWindowShouldClose(engine.window)) {
        engine.setupViewport();
        SagA.draw();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();// ...im Hauptloop, vor ImGui::End()
        static bool bh_type = false;
        const char* types[] = { "Schwarzschild (nicht-rotierend)", "Kerr (rotierend)" };
        
        ImGui::Begin("Controls");
        ImGui::Text("Number of Rays: %d", (int)rays.size());
        ImGui::Text("Black Hole Mass: %.2e kg", SagA.mass);
        ImGui::Text("Schwarzschild Radius: %.2e m", SagA.r_s);
        ImGui::Checkbox("Rotierend blackhole", &bh_type);
        
        if(ImGui::Button("Reset Rays")) {
            rays.clear();
            for (float y = -engine.physHeight; y < engine.physHeight; y += 1e10) {
                rays.push_back(Ray(vec2(-engine.physWidth, y), vec2(c, 0.0)));
            }
        }

        float a = 0.0f; // Kerr-Parameter (a = J/M, Drehimpuls pro Masse)
        if (bh_type == true) {
            //ImGui::SliderFloat("Kerr-Parameter a [m]", &a, 0.0f, SagA.r_s * 0.5f, "%.2e");
            a=SagA.r_s * 0.05f; // Maximale Rotation   
        }
        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
        if (a > SagA.r_s * 0.1f) ImGui::TextColored(ImVec4(1,0.5,0,1), "Warnung: Näherung nur für kleine Kerr-Parameter gültig!");
        ImGui::End();
       for (auto& ray : rays) {
            ray.step(1.0f, SagA.r_s, a, bh_type);
            ray.drawTrailAndPoint(rays);
        }
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(engine.window);

        glfwPollEvents();
    }

    return 0;
}
