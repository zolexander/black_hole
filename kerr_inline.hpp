
#include "blackholesim.hpp"
#include <cmath>
namespace BlackholeSim
{
    inline double Delta(double r, double a)
    {
        return r * r - 2.0 * r + a * a; // geometric units M=1 assumption (simplified)
    }

    // A sample R(r, L, a) placeholder; the real R_of_r is more complex.
    inline double R_of_r(double r, double L, double a)
    {
        // Toy potential: allow radial motion; keep positive in some region.
        double U = 1.0 - 2.0 / r;
        return std::max(0.0, U - (L * L) / (r * r));
    }

    // simplified dr/dlambda and dphi/dlambda functions (toy)
    inline double f_r(double r, double L, double a)
    {
        double val = R_of_r(r, L, a);
        return -sqrt(std::max(0.0, val)) / (r * r + 1e-9);
    }
    inline double f_phi(double r, double L, double a)
    {
        double sig = r * r;
        double del = Delta(r, a);
        double A = (r * r + a * a) - a * L;
        // simplified variant that depends only on r,L,a
        return (a / del) * A - (a - L) / sig;
    }

}