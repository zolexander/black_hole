// Kerr geodesic integration (RK4) implementation
// Keeps a minimal set of includes and preserves the public API.

#include <kerr_inline.hpp>
#include <kerrintegrate.hpp>

namespace BlackholeSim
{

    // Compute time derivatives for the compact state (r, phi)
    // using helpers from kerr_inline.hpp
    static inline KerrState deriv(const KerrState& st, double L, double a_spin)
    {
        return { f_r(st.r, L, a_spin), f_phi(st.r, L, a_spin) };
    }

    KerrState integrateKerr(const KerrState& s0, double L, double& h, double a_spin)
    {
        // Early out for degenerate step sizes
        if (h == 0.0)
        {
            return s0;
        }

        const KerrState k1 = deriv(s0, L, a_spin);
        const KerrState s2 { s0.r + 0.5 * h * k1.r, s0.phi + 0.5 * h * k1.phi };
        const KerrState k2 = deriv(s2, L, a_spin);
        const KerrState s3 { s0.r + 0.5 * h * k2.r, s0.phi + 0.5 * h * k2.phi };
        const KerrState k3 = deriv(s3, L, a_spin);
        const KerrState s4 { s0.r + h * k3.r, s0.phi + h * k3.phi };
        const KerrState k4 = deriv(s4, L, a_spin);

        KerrState out;
        out.r   = s0.r   + (h / 6.0) * (k1.r   + 2.0 * k2.r   + 2.0 * k3.r   + k4.r);
        out.phi = s0.phi + (h / 6.0) * (k1.phi + 2.0 * k2.phi + 2.0 * k3.phi + k4.phi);
        return out;
    }
} // namespace BlackholeSim