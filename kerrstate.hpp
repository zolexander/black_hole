#ifndef BLACKHOLESIM_KERRSTATE_HPP
#define BLACKHOLESIM_KERRSTATE_HPP
// ... dein bisheriger Code ...
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
namespace BlackholeSim
{
    enum class Mode { Kerr = 0, Test = 1 };

    struct KerrState
    {
        double r;
        double phi;
    };
    struct CartesianState
    {
        double x, y;
        double vx, vy;
    };

    struct Photon
    {
        KerrState s;
        double L;
        double h;
        bool alive;
        std::vector<glm::vec2> trail;

        Photon(double L_);

        void reset(double L_, const KerrState &startState);
    };
    struct TestPhoton
    {
        CartesianState s;
        bool alive;
        std::vector<glm::vec2> trail;

        TestPhoton(double x0, double y0);
        void reset(double x0, double y0);
    };
}
#endif