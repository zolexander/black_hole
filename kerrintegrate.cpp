#include <iostream>
#include <blackholesim.hpp>
#include <kerr_inline.hpp>
#include <kerrintegrate.hpp>

namespace BlackholeSim
{

    KerrState integrateKerr(const KerrState &s0, double L, double &h, double a_spin)
    {
        auto f = [&](const KerrState &st) -> KerrState
        {
            return {f_r(st.r, L, a_spin), f_phi(st.r, L, a_spin)};
        };

        KerrState k1 = f(s0);
        KerrState s2 = {s0.r + 0.5 * h * k1.r, s0.phi + 0.5 * h * k1.phi};
        KerrState k2 = f(s2);
        KerrState s3 = {s0.r + 0.5 * h * k2.r, s0.phi + 0.5 * h * k2.phi};
        KerrState k3 = f(s3);
        KerrState s4 = {s0.r + h * k3.r, s0.phi + h * k3.phi};
        KerrState k4 = f(s4);

        KerrState out;
        out.r = s0.r + (h / 6.0) * (k1.r + 2 * k2.r + 2 * k3.r + k4.r);
        out.phi = s0.phi + (h / 6.0) * (k1.phi + 2 * k2.phi + 2 * k3.phi + k4.phi);
        return out;
    }

} // namespace BlackholeSim