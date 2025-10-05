#include "engine.hpp"
#include "blackhole_struct.hpp"
#include "kerr_inline.hpp"
#include "kerrintegrate.hpp"
#include "panel.hpp"
#include "utils.hpp"
#include <GLFW/glfw3.h>
#include <filesystem>
#include <fstream>
#include <optional>
#include <cstring>
namespace BlackholeSim {
Engine::Engine() {}
// ============================================================================
// INITIALIZATION FUNCTIONS
// ============================================================================

/**
 * @brief Initialize OpenGL context, shaders, and ImGui
 * @param title Window title
 * @param w Window width
 * @param h Window height
 * @return true if initialization successful, false otherwise
 */
bool Engine::initGL(const char *title, int w, int h) {
  if (!glfwInit()) {
    std::cerr << "GLFW init failed\n";
    return false;
  }

  // Compute shaders + SSBOs require OpenGL 4.3+
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  window = glfwCreateWindow(w, h, title, nullptr, nullptr);
  if (!window) {
    std::cerr << "glfwCreateWindow failed\n";
    glfwTerminate();
    return false;
  }

  glfwMakeContextCurrent(window);
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
    std::cerr << "gladLoadGLLoader failed\n";
    return false;
  }

  // ImGui init
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO &io = ImGui::GetIO();
  (void)io;
  ImGui::StyleColorsDark();
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init("#version 330");
  shader1 = BlackholeSim::Utils::Shader("photon.vert", "photon.frag");
  shader1.use();
  prog = shader1.ID;
  zoomLoc = shader1.getUniformLocation("zoom");
  projLoc = shader1.getUniformLocation("proj");
  viewLoc = shader1.getUniformLocation("view");

  // VAO/VBO setup
  glGenVertexArrays(1, &vao);
  glGenBuffers(1, &vbo);
  glBindVertexArray(vao);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  // Interleaved layout: [vec2 pos][vec3 color] = 5 floats
  glEnableVertexAttribArray(0); // pos
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void *)0);
  glEnableVertexAttribArray(1); // color
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float),
                        (void *)(2 * sizeof(float)));
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);

  // blending
  glEnable(GL_BLEND);
  // Use standard alpha blending for clearer lines
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDisable(GL_DEPTH_TEST);

  // projection / view
  proj = glm::ortho(-40.0f, 40.0f, -30.0f, 30.0f, -1.0f, 1.0f);
  view = glm::mat4(1.0f);

  // Check GL version supports compute (>= 4.3)
  GLint glMajor = 0, glMinor = 0;
  glGetIntegerv(GL_MAJOR_VERSION, &glMajor);
  glGetIntegerv(GL_MINOR_VERSION, &glMinor);
  bool supportsCompute = (glMajor > 4) || (glMajor == 4 && glMinor >= 3);

  // Initialize photon simulation
  resetPhotons(photonCount);

  // Initialize compute shader only if requested and supported
  if (useGPUPaths && supportsCompute) {
    photonCompute = BlackholeSim::Utils::Shader("compute_photons.comp");
    initPhotonSSBOs();
    uploadPhotonSSBOs();
  } else if (!supportsCompute && useGPUPaths) {
    // Requested GPU but not supported; fall back
    useGPUPaths = false;
    std::cerr << "GPU compute disabled: OpenGL " << glMajor << "." << glMinor
              << " < 4.3. Falling back to CPU path.\n";
  }

  return true;
}

// ============================================================================
// PHOTON MANAGEMENT
// ============================================================================

/**
 * @brief Reset and initialize photons for simulation
 * @param count Number of photons to create
 */
void Engine::resetPhotons(int count) {
  photonCount = count;

  if (mode == Mode::Kerr) {
    initializeKerrPhotons();
  } else {
    initializeTestPhotons();
  }
  if (useGPUPaths) {
    uploadPhotonSSBOs();
  }
}

/**
 * @brief Initialize Kerr photons with proper initial conditions
 */
void Engine::initializeKerrPhotons() {
  kerrPhotons.clear();
  kerrPhotons.reserve(photonCount);

  // Create photons spread along vertical line at initial x position
  const double x0 = Constants::INITIAL_X;
  const double yMin = Constants::Y_MIN;
  const double yMax = Constants::Y_MAX;

  for (int i = 0; i < photonCount; ++i) {
    // Calculate y position with even distribution
    const double y0 = calculatePhotonYPosition(i, yMin, yMax);

    // Convert Cartesian to polar coordinates (r, phi)
    const double r0 = std::sqrt(x0 * x0 + y0 * y0);
    const double phi0 = std::atan2(y0, x0);
    const double L = y0; // Use y-coordinate as impact parameter L

    Photon photon(L);
    photon.reset(L, {r0, phi0});
    photon.trail.reserve(trailLength);
    kerrPhotons.emplace_back(std::move(photon));
  }
}

/**
 * @brief Initialize test photons for Newtonian simulation
 */
void Engine::initializeTestPhotons() {
  testPhotons.clear();
  testPhotons.reserve(photonCount);

  const double x0 = Constants::INITIAL_X;
  const double yMin = Constants::Y_MIN;
  const double yMax = Constants::Y_MAX;

  for (int i = 0; i < photonCount; ++i) {
    const double y0 = calculatePhotonYPosition(i, yMin, yMax);

    TestPhoton photon(x0, y0);
    photon.reset(x0, y0);
    photon.trail.reserve(trailLength);
    testPhotons.emplace_back(std::move(photon));
  }
}

/**
 * @brief Calculate Doppler shift color based on relative velocity
 * @param s0 Initial state
 * @param s1 Final state
 * @return Color vector
 */
glm::vec3 Engine::getDopplerColor(const KerrState &s0, const KerrState &s1) {
  double x0 = s0.r * cos(s0.phi);
  double y0 = s0.r * sin(s0.phi);
  double x1 = s1.r * cos(s1.phi);
  double y1 = s1.r * sin(s1.phi);

  glm::vec2 vel(x1 - x0, y1 - y0);
  if (glm::length(vel) < 1e-8)
    return glm::vec3(1.0f); // neutral (weiß)

  vel = glm::normalize(vel);

  // Beobachterrichtung: +x
  glm::vec2 obsDir(1.0, 0.0);

  double doppler = glm::dot(vel, obsDir);

  glm::vec3 red(1.0, 0.2, 0.2);
  glm::vec3 blue(0.2, 0.2, 1.0);
  glm::vec3 white(1.0, 1.0, 1.0);

  if (doppler >= 0)
    return glm::mix(white, blue, float(doppler));
  else
    return glm::mix(red, white, float(doppler + 1.0));
}
/**
 * @brief Calculate evenly distributed y-position for photon initialization
 * @param index Photon index
 * @param yMin Minimum y value
 * @param yMax Maximum y value
 * @return Calculated y position
 */
double Engine::calculatePhotonYPosition(int index, double yMin,
                                        double yMax) const {
  const int denominator = std::max(1, photonCount - 1);
  return yMin + (yMax - yMin) * (index / static_cast<double>(denominator));
}
// ============================================================================
// PHYSICS UPDATE FUNCTIONS
// ============================================================================

/**
 * @brief Update Kerr photons using relativistic integration
 * @param a_spin Black hole spin parameter
 */
void Engine::updateKerrPhotons(double a_spin) {
  if (useGPUPaths) {
    // When GPU is active, skip CPU integration. Dispatch handled in run().
    return;
  }
  for (auto &photon : kerrPhotons) {
    if (!photon.alive) {
      continue;
    }
    KerrState old = photon.s;
    // Integrate photon trajectory in Kerr coordinates
    photon.s = integrateKerr(photon.s, photon.L, photon.h, a_spin);

    // Convert polar coordinates to Cartesian for rendering
    const float x = static_cast<float>(photon.s.r * std::cos(photon.s.phi));
    const float y = static_cast<float>(photon.s.r * std::sin(photon.s.phi));
    glm::vec3 col = getDopplerColor(old, photon.s);

    photon.trail.emplace_back(x, y, col);

    // Maintain trail length by removing old points
    trimPhotonTrail(photon.trail);

    // Check if photon has crossed the event horizon
    if (photon.s.r < Constants::SCHWARZSCHILD_RADIUS) {
      photon.alive = false;
    }
  }
}

/**
 * @brief Update test photons using Newtonian gravity approximation
 * @param h Integration time step
 * @param M Mass of the central object
 */
void Engine::updateTestPhotons(double h, double M) {
  for (auto &photon : testPhotons) {
    if (!photon.alive) {
      continue;
    }

    // Calculate distance from center
    const double r =
        std::sqrt(photon.s.x * photon.s.x + photon.s.y * photon.s.y);
    const double schwarzschildRadius = Constants::SCHWARZSCHILD_RADIUS * M;

    // Calculate Newtonian gravitational acceleration
    const double rCubed =
        r * r * r + Constants::EPSILON; // Prevent division by zero
    const double ax = -M * photon.s.x / rCubed;
    const double ay = -M * photon.s.y / rCubed;

    // Update velocity
    photon.s.vx += ax * h;
    photon.s.vy += ay * h;

    // Normalize velocity to light speed (c = 1 in our units)
    const double velocityMagnitude =
        std::sqrt(photon.s.vx * photon.s.vx + photon.s.vy * photon.s.vy);
    if (velocityMagnitude > Constants::VELOCITY_THRESHOLD) {
      photon.s.vx /= velocityMagnitude;
      photon.s.vy /= velocityMagnitude;
    }

    // Update position with scaled step for better visibility
    photon.s.x += photon.s.vx * h * Constants::VELOCITY_SCALE;
    photon.s.y += photon.s.vy * h * Constants::VELOCITY_SCALE;

    // Add point to trail
    photon.trail.emplace_back(static_cast<float>(photon.s.x),
                              static_cast<float>(photon.s.y));

    // Maintain trail length
    trimPhotonTrail(photon.trail);

    // Check if photon has been absorbed
    if (r < schwarzschildRadius) {
      photon.alive = false;
    }
  }
}

void Engine::drawKerrPhotons(const std::vector<Photon> &photons, GLuint shader,
                             GLuint vao, GLuint vbo, const glm::mat4 &proj,
                             const glm::mat4 &view, float zoom) {
  shader1.use();
  shader1.setMat4("proj", proj);
  shader1.setMat4("view", view);
  shader1.setFloat("zoom", static_cast<float>(zoom));

  glBindVertexArray(vao);
  for (auto &p : photons) {
    if (p.trail.empty())
      continue;

    std::vector<float> verts;
    verts.reserve(p.trail.size() * 5); // 2 pos + 3 color

    for (const auto &tp : p.trail) {
      verts.push_back(tp.pos.x);
      verts.push_back(tp.pos.y);
      verts.push_back(tp.color.r);
      verts.push_back(tp.color.g);
      verts.push_back(tp.color.b);
    }

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(verts.size() * sizeof(float)),
                 verts.data(), GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float),
                          (void *)0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float),
                          (void *)(2 * sizeof(float)));
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);

    if (p.trail.size() >= 2) {
      glDrawArrays(GL_LINE_STRIP, 0, static_cast<GLsizei>(p.trail.size()));
    } else {
      glPointSize(5.0f);
      glDrawArrays(GL_POINTS, 0, 1);
    }

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
  }

  glBindVertexArray(0);
}
void Engine::drawTestPhotons(GLuint shader, const glm::mat4 &projMat,
                             const glm::mat4 &viewMat) {
  glUseProgram(shader);
  glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(projMat));
  glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(viewMat));
  glUniform1f(zoomLoc, (float)zoom);

  glBindVertexArray(vao);

  for (const auto &p : testPhotons) {
    if (p.trail.size() < 2)
      continue;

    // Build interleaved [pos(2), color(3)] like Kerr path
    std::vector<float> verts;
    verts.reserve(p.trail.size() * 5);

    int vertsCount = (int)p.trail.size();
    const int segs = 3;
    for (int s = 0; s < segs; ++s) {
      int start = (vertsCount * s) / segs;
      int end = (vertsCount * (s + 1)) / segs;
      int count = end - start;
      if (count < 2)
        continue;

      // Choose brightness via alpha weighting (shader has no alpha)
      float alpha = (s == 0 ? Constants::ALPHA_FAINT
                            : (s == 1 ? Constants::ALPHA_MEDIUM
                                      : Constants::ALPHA_BRIGHT));
      const float r = Constants::TEST_RED * alpha;
      const float g = Constants::TEST_COLOR_GREEN * alpha;
      const float b = Constants::TEST_BLUE * alpha;

      // Append this segment's vertices with fixed color
      for (int i = start; i < end; ++i) {
        const glm::vec2 &pt = p.trail[i];
        verts.push_back(pt.x);
        verts.push_back(pt.y);
        verts.push_back(r);
        verts.push_back(g);
        verts.push_back(b);
      }

      glBindBuffer(GL_ARRAY_BUFFER, vbo);
      glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(float), verts.data(),
                   GL_DYNAMIC_DRAW);
      glDrawArrays(GL_LINE_STRIP, 0, count);
      verts.clear();
    }
  }

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
  glUseProgram(0);
}

void Engine::run() {
  if (!window)
    return;

  double lastTime = glfwGetTime();
  while (!glfwWindowShouldClose(window)) {
    glfwPollEvents();
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    Controlpanel(*this);

    // Update physics simulation
    double now = glfwGetTime();
    float dt = static_cast<float>(now - lastTime);
    lastTime = now;
    if (mode == Mode::Kerr) {
      if (useGPUPaths) {
        // speed_scale and doppler_strength are artistic controls
        // Moderate advance per frame for stable visuals
        dispatchPhotonCompute(dt > 0 ? dt : 0.016f, static_cast<float>(a_spin), 10.0f, 0.5f);
        readTrailsBackToCPU();
      } else {
        updateKerrPhotons(a_spin);
      }
    } else {
      updateTestPhotons(Constants::INTEGRATION_STEP, Constants::MASS_DEFAULT);
    }

    // Periodic debug output for monitoring simulation state
    outputDebugInfo();

    // render
    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClearColor(0.02f, 0.02f, 0.03f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    // draw photons
    if (mode == Mode::Kerr) {
      drawKerrPhotons(kerrPhotons, prog, vao, vbo, proj, view, zoom);
      double rp = bh.r_plus();
      double rm = bh.r_minus();
      if (std::isfinite(rp)) {
        BlackholeSim::drawBlackHoleVisuals(bh, shader1, vao, vbo, proj, view,
                                           (float)zoom, showHorizons,
                                           showErgosphere);
      }
      if (std::isfinite(rm)) {
        BlackholeSim::drawBlackHoleVisuals(bh, shader1, vao, vbo, proj, view,
                                           (float)zoom, showHorizons,
                                           showErgosphere);
      }
    } else {
      // Test mode rendering
      drawTestPhotons(prog, proj, view);
    }
    // Finalize Dear ImGui for this frame (builds the command lists)
    ImGui::Render();

    // Reset GL state to what ImGui expects to avoid crashes in the backend
    glUseProgram(0);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendEquation(GL_FUNC_ADD);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    glfwSwapBuffers(window);
  }
}

/**
 * @brief Output periodic debug information about simulation state
 */
void Engine::outputDebugInfo() {
  static int debugCounter = 0;
  if ((debugCounter++ % Constants::DEBUG_OUTPUT_FREQUENCY) == 0) {
    if (!kerrPhotons.empty()) {
      const auto &photon = kerrPhotons.front();
      std::cerr << "DEBUG[Kerr] L=" << photon.L << " r=" << photon.s.r
                << " phi=" << photon.s.phi
                << " trail_size=" << photon.trail.size()
                << " alive=" << (photon.alive ? "true" : "false") << std::endl;
    }
    if (!testPhotons.empty()) {
      const auto &photon = testPhotons.front();
      std::cerr << "DEBUG[Test] x=" << photon.s.x << " y=" << photon.s.y
                << " trail_size=" << photon.trail.size()
                << " alive=" << (photon.alive ? "true" : "false") << std::endl;
    }
  }
}

// ----------------------------------------------------------------------------
// GPU compute helpers (definitions)
// ----------------------------------------------------------------------------
void Engine::initPhotonSSBOs() {
  if (ssboPhotons == 0) photonCompute.GenBuffers(1, &ssboPhotons);
  if (ssboTrails == 0) photonCompute.GenBuffers(1, &ssboTrails);

  // Allocate initial buffers
  GLsizeiptr photonsSize = static_cast<GLsizeiptr>(photonCount) * 48; // std430 struct size
  photonCompute.BindBuffer(GL_SHADER_STORAGE_BUFFER, ssboPhotons);
  photonCompute.BufferData(GL_SHADER_STORAGE_BUFFER, photonsSize, nullptr, GL_DYNAMIC_DRAW);
  photonCompute.bindSSBO(0, ssboPhotons);

  GLsizeiptr trailsSize = static_cast<GLsizeiptr>(photonCount) * gpuTrailLen * sizeof(glm::vec2);
  
  photonCompute.BindBuffer(GL_SHADER_STORAGE_BUFFER, ssboTrails);
  photonCompute.BufferData(GL_SHADER_STORAGE_BUFFER, trailsSize, nullptr, GL_DYNAMIC_DRAW);
  photonCompute.bindSSBO(1, ssboTrails);
  photonCompute.BindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
}

void Engine::uploadPhotonSSBOs() {
  if (ssboPhotons == 0 || ssboTrails == 0) {
    initPhotonSSBOs();
  }
  struct PhotonGPU {
    float s_x, s_y; // r, phi
    float L;
    float h;
    float alive;
    float pad_[3];
    float col[4];
  };

  std::vector<PhotonGPU> data;
  data.resize(photonCount);
  for (int i = 0; i < photonCount; ++i) {
    PhotonGPU p{};
    if (i < static_cast<int>(kerrPhotons.size())) {
      const auto &src = kerrPhotons[i];
      p.s_x = static_cast<float>(src.s.r);
      p.s_y = static_cast<float>(src.s.phi);
      p.L = static_cast<float>(src.L);
      p.h = static_cast<float>(src.h > 0 ? src.h : 0.02);
      p.alive = src.alive ? 1.0f : 0.0f;
      p.col[0] = 1.0f; p.col[1] = 1.0f; p.col[2] = 1.0f; p.col[3] = 1.0f;
    } else {
      p.alive = 0.0f;
    }
    data[i] = p;
  }

  // Upload photons
  GLsizeiptr photonsSize = static_cast<GLsizeiptr>(data.size() * sizeof(PhotonGPU));
  photonCompute.BindBuffer(GL_SHADER_STORAGE_BUFFER, ssboPhotons);
  photonCompute.BufferData(GL_SHADER_STORAGE_BUFFER, photonsSize, data.data(), GL_DYNAMIC_DRAW);
  photonCompute.bindSSBO(0, ssboPhotons);

  // Ensure trails buffer sized
  GLsizeiptr trailsSize = static_cast<GLsizeiptr>(photonCount) * gpuTrailLen * sizeof(glm::vec2);
  photonCompute.BindBuffer(GL_SHADER_STORAGE_BUFFER, ssboTrails);
  // Pre-fill trails with current states so we don't draw zeros
  std::vector<glm::vec2> initTrails;
  initTrails.resize(static_cast<size_t>(photonCount) * static_cast<size_t>(gpuTrailLen));
  for (int i = 0; i < photonCount; ++i) {
    float r = (i < (int)kerrPhotons.size()) ? static_cast<float>(kerrPhotons[i].s.r) : 0.0f;
    float phi = (i < (int)kerrPhotons.size()) ? static_cast<float>(kerrPhotons[i].s.phi) : 0.0f;
    glm::vec2 s(r, phi);
    for (int t = 0; t < gpuTrailLen; ++t) {
      initTrails[i * gpuTrailLen + t] = s;
    }
  }
  photonCompute.BufferData(GL_SHADER_STORAGE_BUFFER, trailsSize, initTrails.data(), GL_DYNAMIC_DRAW);
  photonCompute.bindSSBO(1, ssboTrails);
  photonCompute.BindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
}

void Engine::dispatchPhotonCompute(float dt, float aSpin, float speedScale, float dopplerStrength) {
  photonCompute.use();
  photonCompute.bindSSBO(0,ssboPhotons);
  photonCompute.bindSSBO(1, ssboTrails);
  // Uniforms matching compute_photons.comp
  photonCompute.setFloat("a_spin", aSpin);
  photonCompute.setFloat("speed_scale", speedScale);
  photonCompute.setFloat("doppler_strength", dopplerStrength);
  photonCompute.setInt("photon_count", photonCount);
  photonCompute.setInt("trailLen", gpuTrailLen);
  photonCompute.setInt("frameIndex", frameIndex);

  GLuint groups = static_cast<GLuint>((photonCount + 255) / 256);
  photonCompute.DispatchCompute(groups,1,1);
  frameIndex++;
}

void Engine::readTrailsBackToCPU() {
  if (kerrPhotons.size() < static_cast<size_t>(photonCount)) return;

  struct PhotonGPU {
    float s_x, s_y; // r, phi
    float L;
    float h;
    float alive;
    float pad_[3];
    float col[4];
  };

  // Read photons
  
  glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboPhotons);
  PhotonGPU *pData = (PhotonGPU*)glMapBufferRange(GL_SHADER_STORAGE_BUFFER, 0,
                         photonCount * sizeof(PhotonGPU), GL_MAP_READ_BIT);
  std::vector<PhotonGPU> photonsCPU;
  photonsCPU.resize(photonCount);
  if (pData) {
    std::memcpy(photonsCPU.data(), pData, photonCount * sizeof(PhotonGPU));
    glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
  }

  // Read trails
  glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboTrails);
  const size_t trailsCount = static_cast<size_t>(photonCount) * static_cast<size_t>(gpuTrailLen);
  glm::vec2 *tData = (glm::vec2*)glMapBufferRange(GL_SHADER_STORAGE_BUFFER, 0,
                       trailsCount * sizeof(glm::vec2), GL_MAP_READ_BIT);
  if (!tData) {
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    return;
  }

  for (int i = 0; i < photonCount; ++i) {
    auto &cpuPhoton = kerrPhotons[i];
    const PhotonGPU &pg = photonsCPU[i];
    glm::vec3 col(pg.col[0], pg.col[1], pg.col[2]);
    // Append a few of the most recent ring entries to grow trail faster
    if (frameIndex > 0u) {
      unsigned maxAppend = 5u;
      unsigned available = std::min<unsigned>(static_cast<unsigned>(gpuTrailLen), frameIndex);
      unsigned count = std::min<unsigned>(maxAppend, available);
      for (unsigned j = 0; j < count; ++j) {
        unsigned idx = (frameIndex - 1u - j) % static_cast<unsigned>(gpuTrailLen);
        const glm::vec2 s = tData[i * gpuTrailLen + static_cast<int>(idx)]; // r, phi
        float r = s.x;
        float phi = s.y;
        if ((r > 0.0f) && std::isfinite(r) && std::isfinite(phi)) {
          float x = r * std::cos(phi);
          float y = r * std::sin(phi);
          cpuPhoton.trail.emplace_back(x, y, col);
        }
      }
    }
    // Trim to configured trail length
    if (cpuPhoton.trail.size() > trailLength) {
      const size_t removeCount = cpuPhoton.trail.size() - trailLength;
      cpuPhoton.trail.erase(cpuPhoton.trail.begin(), cpuPhoton.trail.begin() + removeCount);
    }
    // Update state for debug
    cpuPhoton.s.r = pg.s_x;
    cpuPhoton.s.phi = pg.s_y;
    cpuPhoton.h = pg.h;
    cpuPhoton.L = pg.L;
    cpuPhoton.alive = (pg.alive > 0.5f);
  }

  glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
  glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
}

Engine::~Engine() {
  if (vbo)
    glDeleteBuffers(1, &vbo);
  if (vao)
    glDeleteVertexArrays(1, &vao);
  if (prog)
    glDeleteProgram(prog);
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();
}
}