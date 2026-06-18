//
// Created by ozdalkiran-l on 1/8/26.
//

#include "GaugeField.h"

#include "../su3/utils.h"

// Initialises the gauge field with random SU3 matrices
void GaugeField::hot_start(const Geometry& geo, std::mt19937_64& rng) {
    for (int t = 0; t < T; t++)
        for (int z = 1; z <= L_int; z++)
            for (int y = 1; y <= L_int; y++)
                for (int x = 1; x <= L_int; x++) {
                    size_t site = geo.index(x, y, z, t);
                    for (int mu = 0; mu < 4; mu++) {
                        if (t == T - 1 and mu == 3)
                            view_link(site, mu) = SU3::Identity();
                        else
                            view_link(site, mu) = random_su3(rng);
                    }
                }
}

// Initialises the gauge field with identity matrices
void GaugeField::cold_start(const Geometry& geo) {
    for (int t = 0; t < T; t++)
        for (int z = 1; z <= L_int; z++)
            for (int y = 1; y <= L_int; y++)
                for (int x = 1; x <= L_int; x++) {
                    size_t site = geo.index(x, y, z, t);
                    for (int mu = 0; mu < 4; mu++) {
                        view_link(site, mu) = SU3::Identity();
                    }
                }
}

// Projects a link on SU3 using Gramm-Schmidt
void GaugeField::projection_su3(size_t site, int mu) {
    auto U = view_link(site, mu);

    SU3 temp = U;

    Eigen::Vector3cd c0 = temp.col(0);
    c0.normalize();

    Eigen::Vector3cd c1 = temp.col(1);
    c1 -= c0 * c0.dot(c1);
    c1.normalize();

    Eigen::Vector3cd c2;
    c2(0) = std::conj(c0(1) * c1(2) - c0(2) * c1(1));
    c2(1) = std::conj(c0(2) * c1(0) - c0(0) * c1(2));
    c2(2) = std::conj(c0(0) * c1(1) - c0(1) * c1(0));

    temp.col(0) = c0;
    temp.col(1) = c1;
    temp.col(2) = c2;

    U = temp;
}

// Projects the whole field on SU3
void GaugeField::project_field_su3(const Geometry& geo) {
    for (int t = 0; t < T; t++)
        for (int z = 1; z <= L_int; z++)
            for (int y = 1; y <= L_int; y++)
                for (int x = 1; x <= L_int; x++) {
                    size_t site = geo.index(x, y, z, t);
                    for (int mu = 0; mu < 4; mu++) {
                        projection_su3(site, mu);
                    }
                }
}

// Computes the sum of all staples of a site
void GaugeField::compute_staple(const Geometry& geo, size_t site, int mu, SU3& staple) const {
    staple.setZero();
    size_t link_idx = site * 4 + mu;

    // Forward staples
    for (size_t i = geo.fwd_start[link_idx]; i < geo.fwd_start[link_idx + 1]; ++i) {
        const auto& s = geo.fwd_staples_opt[i];
        staple += s.coeff * (view_link_const_off(s.off0) * 
                            view_link_const_off(s.off1).adjoint() * 
                            view_link_const_off(s.off2).adjoint());
    }

    // Backward staples
    for (size_t i = geo.bwd_start[link_idx]; i < geo.bwd_start[link_idx + 1]; ++i) {
        const auto& s = geo.bwd_staples_opt[i];
        staple += s.coeff * (view_link_const_off(s.off0).adjoint() * 
                            view_link_const_off(s.off1).adjoint() * 
                            view_link_const_off(s.off2));
    }
}

// Performs a random gauge_transform on the field
void GaugeField::gauge_transform(const Geometry& geo, std::mt19937_64& rng) {
    GaugeField transform(geo);
    transform.hot_start(geo, rng);
    for (int t = 0; t < T; t++)
        for (int z = 1; z <= L_int; z++)
            for (int y = 1; y <= L_int; y++)
                for (int x = 1; x <= L_int; x++) {
                    size_t site = geo.index(x, y, z, t);
                    for (int mu = 0; mu < 4; mu++) {
                        if (t < T - 1 or mu < 3) //ie if not a OBC fixed link
                            view_link(site, mu) =
                                transform.view_link_const(site, mu) * view_link_const(site, mu) *
                                transform.view_link_const(geo.get_neigh(site, mu, up), mu)
                                    .adjoint();
                    }
                }
};
