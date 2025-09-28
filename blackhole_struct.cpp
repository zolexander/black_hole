#include "blackhole_struct.hpp"
#include <cmath>
#include <functional>
#include <glad/glad.h>
#include <glm/glm.hpp>
#include <vector>
#include <iostream>

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
  
    // Set color/alpha and issue the draw call
    glUniform3fv(colorLoc, 1, &color[0]);
    glUniform1f(alphaLoc, alpha);
    glDrawArrays(GL_LINE_STRIP, 0, static_cast<GLsizei>(points.size()));
  
    // Restore previous caller state
    glBindBuffer(GL_ARRAY_BUFFER, prevArrayBuffer);
    glBindVertexArray(prevVAO);
    BH_DBG("[DBG] drawPolylineFromVec: end");
  }
// -------------------------------------------------------------------------------------
// High-level: draw black hole horizons and ergosphere as colored polylines
// Notes:
//  - Does minimal GL state changes (keeps global state for ImGui stability).
//  - Expects shader program to provide uniforms: proj, view, zoom, color, alpha.
// -------------------------------------------------------------------------------------
void BlackholeSim::drawBlackHoleVisuals(
    const BlackHole &bh,
    GLuint shader,
    GLuint vao,
    GLuint vbo,
    const glm::mat4 &proj,
    const glm::mat4 &view,
    float zoom,
    bool showHorizons,
    bool showErgosphere)
  {
    BH_DBG("[DBG] drawBlackHoleVisuals: begin");

    // Save/restore only the active program to avoid interfering with ImGui
    GLint prevProgram = 0;
    glGetIntegerv(GL_CURRENT_PROGRAM, &prevProgram);
    glUseProgram(shader);

    if (!glIsProgram(shader)) { BH_DBG("GL ERROR: Invalid shader program in drawBlackHoleVisuals"); return; }

    // Uniforms
    const GLint projLoc  = glGetUniformLocation(shader, "proj");
    const GLint viewLoc  = glGetUniformLocation(shader, "view");
    const GLint zoomLoc  = glGetUniformLocation(shader, "zoom");
    const GLint colorLoc = glGetUniformLocation(shader, "color");
    const GLint alphaLoc = glGetUniformLocation(shader, "alpha");

    glUniformMatrix4fv(projLoc, 1, GL_FALSE, &proj[0][0]);
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, &view[0][0]);
    glUniform1f(zoomLoc, zoom);

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
    glUseProgram(prevProgram);
  }
