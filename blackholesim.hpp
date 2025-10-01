#ifndef BLACKHOLESIM_HPP
#define BLACKHOLESIM_HPP
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
// namespace BlackholeSim
namespace BlackholeSim
{
    enum class Mode
    {
        Kerr = 0,
        Test = 1
    };

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
    struct TrailPoint {
        glm::vec2 pos;
        glm::vec3 color;
    
        TrailPoint(float x, float y, glm::vec3 c) 
            : pos(x, y), color(c) {}
    };
    struct Photon
    {
        KerrState s;
        double L;
        double h;
        bool alive;
        std::vector<TrailPoint> trail;

        Photon(double L_ = 0.0);

        void reset(double L_, const KerrState &startState);
    };
    struct TestPhoton
    {
        CartesianState s;
        bool alive;
        std::vector<glm::vec2> trail;

        TestPhoton(double x0 = -20.0, double y0 = 0.0);
        void reset(double x0, double y0);
    };
}
#endif