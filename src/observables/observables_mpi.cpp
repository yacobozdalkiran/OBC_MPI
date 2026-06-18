#include "observables_mpi.h"

#include <omp.h>
#include <iostream>

// Computation of global mean plaquette with halos embedded in field (needs field halos exchanges
// first)
double mpi::observables::mean_plaquette_global(GaugeField& field, const Geometry& geo,
                                               MpiTopology& topo, int T_min, int T_max) {
    double local_mean_plaquette = mean_plaquette_weighted(field, geo, T_min, T_max);
    double global_mean_plaquette = 0.0;
    MPI_Allreduce(&local_mean_plaquette, &global_mean_plaquette, 1, MPI_DOUBLE, MPI_SUM,
                  topo.cart_comm);
    global_mean_plaquette /= topo.size;
    return global_mean_plaquette;
}

// Computation of levi-civita tensor
int mpi::observables::levi_civita(int mu, int nu, int rho, int sigma) {
    // A stocker dans un tableau
    if (mu == nu || mu == rho || mu == sigma || nu == rho || nu == sigma || rho == sigma) return 0;
    int inv = 0;
    inv += (mu > nu);
    inv += (mu > rho);
    inv += (mu > sigma);
    inv += (nu > rho);
    inv += (nu > sigma);
    inv += (rho > sigma);

    return (inv & 1) ? -1 : +1;  // Parity test on inv
}

// Computes G_{\mu\nu}(site) clover
SU3 mpi::observables::clover_site(const GaugeField& field, const Geometry& geo, size_t site,
                                  int mu, int nu) {
    if (mu == nu) std::cerr << "mu = nu => G(site, mu, nu) = 0\n";
    size_t x = site;
    size_t xpmu = geo.get_neigh(x, mu, up);          // x+mu
    size_t xpnu = geo.get_neigh(x, nu, up);          // x+nu
    size_t xmmu = geo.get_neigh(x, mu, down);        // x-mu
    size_t xmnu = geo.get_neigh(x, nu, down);        // x-nu
    size_t xmmupnu = geo.get_neigh(xmmu, nu, up);    // x-mu+nu
    size_t xmmumnu = geo.get_neigh(xmmu, nu, down);  // x-mu-nu
    size_t xpmumnu = geo.get_neigh(xpmu, nu, down);  // x+mu-nu
    SU3 clover = SU3::Zero();
    clover += field.view_link_const(x, mu) * field.view_link_const(xpmu, nu) *
              field.view_link_const(xpnu, mu).adjoint() * field.view_link_const(x, nu).adjoint();
    clover += field.view_link_const(x, nu) * field.view_link_const(xmmupnu, mu).adjoint() *
              field.view_link_const(xmmu, nu).adjoint() * field.view_link_const(xmmu, mu);
    clover += field.view_link_const(xmmu, mu).adjoint() *
              field.view_link_const(xmmumnu, nu).adjoint() * field.view_link_const(xmmumnu, mu) *
              field.view_link_const(xmnu, nu);
    clover += field.view_link_const(xmnu, nu).adjoint() * field.view_link_const(xmnu, mu) *
              field.view_link_const(xpmumnu, nu) * field.view_link_const(x, mu).adjoint();
    SU3 F = (clover - clover.adjoint()).eval();
    F *= 0.125;
    return F;
}
