#include "blackhole_struct.hpp"
#include <cmath>
#include <functional>
#include <glad/glad.h>
#include <glm/glm.hpp>
#include <vector>
#include <iostream>
#include "utils.hpp"

// -------------------------------------------------------------------------------------
// Geometry helper: sample a parametric radius r(theta) into a closed polyline
// -------------------------------------------------------------------------------------
std::vector<glm::vec2> BlackholeSim::makeEllipsePoints(std::function<double(double)> r_of_theta, int segments = 180)
{
  std::vector<glm::vec2> pts;
  pts.reserve(static_cast<size_t>(segments) + 1);
  for (int i = 0; i <= segments; ++i) {
    const double theta = 2.0 * M_PI * i / segments;
    double r = r_of_theta(theta);
    if (!std::isfinite(r)) r = 0.0;
    const double x = r * std::cos(theta);
    const double y = r * std::sin(theta);
    pts.emplace_back(static_cast<float>(x), static_cast<float>(y));
  }
  return pts;
}
// -------------------------------------------------------------------------------------
// Rendering helper: draw a 2D polyline from a point list using an already configured VAO
// Requirements:
//  - The provided VAO must have attribute 0 as vec2 position bound to ARRAY_BUFFER.
//  - The provided shader program must be in use and have uniform locations valid.
//  - This function keeps GL state impact minimal and restores previous VAO/VBO bindings.
// -------------------------------------------------------------------------------------
void BlackholeSim::drawPolylineFromVec(
    GLuint vao,
    GLuint vbo,
    const std::vector<glm::vec2> &points,
    GLint colorLoc,
    GLint alphaLoc,
    const glm::vec3 &color,
    float alpha)
  {
    if (points.empty()) return;
    BH_DBG("[DBG] drawPolylineFromVec: begin, points=" << points.size());

    if (!glIsVertexArray(vao)) { BH_DBG("GL ERROR: Invalid VAO in drawPolylineFromVec"); return; }
    if (!glIsBuffer(vbo))      { BH_DBG("GL ERROR: Invalid VBO in drawPolylineFromVec"); return; }

    // Backup previous bindings
    GLint prevVAO = 0;
    GLint prevArrayBuffer = 0;
    glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &prevVAO);
    glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &prevArrayBuffer);

    // Upload vertex data
    glBindVertexArray(vao);            BH_DBG("[DBG] drawPolylineFromVec: bound VAO");
    glBindBuffer(GL_ARRAY_BUFFER, vbo);BH_DBG("[DBG] drawPolylineFromVec: bound VBO");
    glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(glm::vec2), points.data(), GL_DYNAMIC_DRAW);
    
    // Configure attribute 0 for tight-packed vec2 positions (override interleaved layout)
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), (void*)0);

    // Use a constant vertex attribute for color (attribute 1) so current shader's inColor works
    glDisableVertexAttribArray(1);
    glVertexAttrib3f(1, color.r * alpha, color.g * alpha, color.b * alpha);
    
    // Issue the draw call
    glDrawArrays(GL_LINE_STRIP, 0, static_cast<GLsizei>(points.size()));

    // Restore previous caller state
    // Restore attribute array 1 (re-enable for interleaved buffers used elsewhere)
    glEnableVertexAttribArray(1);
    // Restore attribute 0 layout back to interleaved [vec2 pos][vec3 color] = 5 floats
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
    glBindBuffer(GL_ARRAY_BUFFER, prevArrayBuffer);
    glBindVertexArray(prevVAO);
    BH_DBG("[DBG] drawPolylineFromVec: end");
  }
// -------------------------------------------------------------------------------------
// Notes:
//  - Does minimal GL state changes (keeps global state for ImGui stability).
//  - Expects shader program to provide uniforms: proj, view, zoom, color, alpha.
// -------------------------------------------------------------------------------------
void BlackholeSim::drawBlackHoleVisuals(
    const BlackHole &bh,
    BlackholeSim::Utils::Shader &shader,
    GLuint vao,
    GLuint vbo,
    const glm::mat4 &proj,
    const glm::mat4 &view,
    float zoom,
    bool showHorizons,
    bool showErgosphere)
  {
    BH_DBG("[DBG] drawBlackHoleVisuals: begin");

    // Bind program via wrapper and set common uniforms
    shader.use();
    shader.setMat4("proj", proj);
    shader.setMat4("view", view);
    shader.setFloat("zoom", zoom);

    // This shader does not use 'color'/'alpha' uniforms; avoid querying to prevent warnings
    const GLint colorLoc = -1;
    const GLint alphaLoc = -1;

    // Sample/draw shapes based on flags
    if (showHorizons) {
      const auto circEH = makeEllipsePoints([&](double){ return bh.r_plus(); }, 256);
      const auto circIH = makeEllipsePoints([&](double){ return bh.r_minus(); }, 256);
      BH_DBG("Horizons: rp=" << bh.r_plus() << " rm=" << bh.r_minus());
      // Outer horizon (red)
      BH_DBG("[DBG] drawBlackHoleVisuals: draw EH");
      drawPolylineFromVec(vao, vbo, circEH, colorLoc, alphaLoc, glm::vec3(1.0f, 0.0f, 0.0f), 1.0f);

      // Inner horizon (green)
      BH_DBG("[DBG] drawBlackHoleVisuals: draw IH");
      drawPolylineFromVec(vao, vbo, circIH, colorLoc, alphaLoc, glm::vec3(0.0f, 1.0f, 0.0f), 1.0f);
    }

    if (showErgosphere) {
      const auto circErgo = makeEllipsePoints([&](double th){ return bh.r_erg(th); }, 256);
      BH_DBG("Ergo: ergo@pi/2=" << bh.r_erg(M_PI/2));
      // Ergosphere (blue, slightly transparent)
      BH_DBG("[DBG] drawBlackHoleVisuals: draw ERGO");
      drawPolylineFromVec(vao, vbo, circErgo, colorLoc, alphaLoc, glm::vec3(0.0f, 0.0f, 1.0f), 0.7f);
    }

    BH_DBG("[DBG] drawBlackHoleVisuals: end");
  }
