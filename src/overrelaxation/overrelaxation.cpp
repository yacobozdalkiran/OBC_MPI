#include "overrelaxation.h"

#include <omp.h>

// Calculates the overrelaxation step for a single SU(2) subgroup block w
// It computes r = (v^\dagger)^2 = v^\dagger v^\dagger where v = w / |w| is in SU(2).
// This is deterministic and preserves the trace (Re Tr(r w) = Re Tr(w)).
SU2q overrelaxation::su2_step(const SU2q& w) {
    double k_sq = w.squaredNorm();
    if (k_sq < 1e-14) {
        return SU2q(1.0, 0.0, 0.0, 0.0);
    }
    double coeff = 2.0 / k_sq;
    double w0 = w(0);
    return SU2q(coeff * w0 * w0 - 1.0, -coeff * w0 * w(1), -coeff * w0 * w(2), -coeff * w0 * w(3));
}

// Performs an overrelaxation update on a single link (sequential)
void overrelaxation::hit(GaugeField& field, const Geometry& geo, size_t site, int mu, SU3& A) {
    field.compute_staple(geo, site, mu, A);

    // Subgroup (0,1)
    SU3 W = field.view_link_const(site, mu) * A;
    SU2q wq = su2_to_quaternion(W.block<2, 2>(0, 0));
    SU2q r = su2_step(wq);
    SU3 R = su2_quaternion_to_su3(r, 0, 1);
    field.view_link(site, mu) = R * field.view_link(site, mu);

    // Subgroup (1,2)
    W = field.view_link_const(site, mu) * A;
    wq = su2_to_quaternion(W.block<2, 2>(1, 1));
    r = su2_step(wq);
    R = su2_quaternion_to_su3(r, 1, 2);
    field.view_link(site, mu) = R * field.view_link(site, mu);

    // Subgroup (0,2)
    W = field.view_link_const(site, mu) * A;
    SU2 wsu2;
    wsu2 << W(0, 0), W(0, 2), W(2, 0), W(2, 2);
    wq = su2_to_quaternion(wsu2);
    r = su2_step(wq);
    R = su2_quaternion_to_su3(r, 0, 2);
    field.view_link(site, mu) = R * field.view_link(site, mu);
}

// Full overrelaxation sweep (sequential/single-thread)
void overrelaxation::sweep(GaugeField& field, const Geometry& geo) {
    SU3 A;
    for (int t = 0; t < geo.T; t++) {
        for (int z = 1; z <= geo.L_int; z++) {
            for (int y = 1; y <= geo.L_int; y++) {
                for (int x = 1; x <= geo.L_int; x++) {
                    size_t site = geo.index(x, y, z, t);
                    for (int mu = 0; mu < 4; mu++) {
                        if (t == geo.T - 1 and mu == 3) continue;
                        hit(field, geo, site, mu, A);
                    }
                }
            }
        }
    }
}

// Performs an overrelaxation hit on a single link (used by checkerboard)
void mpi::overrelaxationcb::hit(GaugeField& field, const Geometry& geo, size_t site, int mu,
                                SU3& A) {
    field.compute_staple(geo, site, mu, A);

    // Subgroup (0,1)
    SU3 W = field.view_link_const(site, mu) * A;
    SU2q wq = su2_to_quaternion(W.block<2, 2>(0, 0));
    SU2q r = ::overrelaxation::su2_step(wq);
    SU3 R = su2_quaternion_to_su3(r, 0, 1);
    field.view_link(site, mu) = R * field.view_link(site, mu);

    // Subgroup (1,2)
    W = field.view_link_const(site, mu) * A;
    wq = su2_to_quaternion(W.block<2, 2>(1, 1));
    r = ::overrelaxation::su2_step(wq);
    R = su2_quaternion_to_su3(r, 1, 2);
    field.view_link(site, mu) = R * field.view_link(site, mu);

    // Subgroup (0,2)
    W = field.view_link_const(site, mu) * A;
    SU2 wsu2;
    wsu2 << W(0, 0), W(0, 2), W(2, 0), W(2, 2);
    wq = su2_to_quaternion(wsu2);
    r = ::overrelaxation::su2_step(wq);
    R = su2_quaternion_to_su3(r, 0, 2);
    field.view_link(site, mu) = R * field.view_link(site, mu);
}

// Performs an overrelaxation sweep on links of the specified parity
void mpi::overrelaxationcb::sweep(GaugeField& field, const Geometry& geo,
                                  site_parity update_parity) {
    for (int mu = 0; mu < 4; mu++) {
#pragma omp parallel for collapse(4)
        for (int t = 0; t < geo.T; t++) {
            for (int z = 1; z <= geo.L_int; z++) {
                for (int y = 1; y <= geo.L_int; y++) {
                    for (int x = 1; x <= geo.L_int; x++) {
                        size_t site = geo.index(x, y, z, t);
                        if (!geo.is_frozen(site, mu) and (geo.get_parity(site) == update_parity) and
                            not(t == geo.T - 1 and mu == 3)) {
                            SU3 A;
                            hit(field, geo, site, mu, A);
                        }
                    }
                }
            }
        }
    }
}

// Runs a specified number of overrelaxation sweeps using checkerboard updates
void mpi::overrelaxationcb::sample(GaugeField& field, const Geometry& geo, int N_sweeps) {
    for (int s = 0; s < N_sweeps; s++) {
        sweep(field, geo, site_parity::EV);
        sweep(field, geo, site_parity::OD);
    }
}
