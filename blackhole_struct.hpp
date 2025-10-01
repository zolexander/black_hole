#pragma once
// =====================================================================================
// Black hole helpers and simple rendering utilities
// =====================================================================================

#ifndef BLACKHOLE_STRUCT_H
#define BLACKHOLE_STRUCT_H

#include <cmath>
#include <functional>
#include <glad/glad.h>
#include <glm/glm.hpp>
#include <vector>
#include "utils.hpp"

// Toggle debug prints in this header (0 = off, 1 = on)
#ifndef BH_STRUCT_DEBUG
#define BH_STRUCT_DEBUG 0
#endif

#if BH_STRUCT_DEBUG
#define BH_DBG(msg) (std::cerr << msg << '\n')
#else
#define BH_DBG(msg) ((void)0)
#endif

namespace BlackholeSim {

// -------------------------------------------------------------------------------------
// Model: Simple Kerr black hole parameters and characteristic radii
// -------------------------------------------------------------------------------------
struct BlackHole {
  // Mass (geometric units) and spin parameter (|a| <= M)
  double M = 1.0;
  double a = 0.7;

  // Outer horizon radius r_+ = M + sqrt(M^2 - a^2) (NaN for naked singularity)
  double r_plus() const {
    const double disc = M * M - a * a;
    if (disc < 0.0)
      return NAN;
    return M + std::sqrt(disc);
  }

  // Inner horizon radius r_- = M - sqrt(M^2 - a^2) (NaN for naked singularity)
  double r_minus() const {
    const double disc = M * M - a * a;
    if (disc < 0.0)
      return NAN;
    return M - std::sqrt(disc);
  }

  // Ergosphere radius (equatorial is r_erg(pi/2) = M + sqrt(M^2 - a^2))
  double r_erg(double theta) const {
    const double c = std::cos(theta);
    const double disc = M * M - a * a * c * c;
    if (disc < 0.0)
      return NAN;
    return M + std::sqrt(disc);
  }
};

std::vector<glm::vec2>
makeEllipsePoints(std::function<double(double)> r_of_theta, int segments);

void drawPolylineFromVec(GLuint vao, GLuint vbo,
                         const std::vector<glm::vec2> &points, GLint colorLoc,
                         GLint alphaLoc, const glm::vec3 &color, float alpha);
void drawBlackHoleVisuals(const BlackHole &bh, BlackholeSim::Utils::Shader &shader, GLuint vao,
                          GLuint vbo, const glm::mat4 &proj,
                          const glm::mat4 &view, float zoom, bool showHorizons,
                          bool showErgosphere);
} // namespace BlackholeSim
#endif