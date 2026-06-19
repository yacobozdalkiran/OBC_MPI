#ifndef INC_4D_MPI_OBSERVABLES_MPI_H
#define INC_4D_MPI_OBSERVABLES_MPI_H

#include "../gauge/GaugeField.h"
#include "../mpi/MpiTopology.h"
#include "../flow/gradient_flow.h"

namespace mpi::observables {
double mean_plaquette_global(GaugeField& field, const Geometry& geo, MpiTopology& topo, int T_min,
                             int T_max);
SU3 clover_site(const GaugeField& field, const Geometry& geo, size_t site, int mu, int nu);
int levi_civita(int mu, int nu, int rho, int sigma);
std::pair<double, double> local_q_e_clover(const GaugeField& field, const Geometry& geo,
                                           size_t site);
std::pair<double, double> topo_q_e_clover(const GaugeField& field, const Geometry& geo);
std::pair<double, double> topo_q_e_clover_global(const GaugeField& field, const Geometry& geo,
                                                 MpiTopology& topo);
std::vector<double> topo_charge_flowed(GaugeField& field, const Geometry& geo, GradientFlow& gf, mpi::MpiTopology& topo, int N_steps_gf, int N_rk_steps);
}  // namespace mpi::observables

#endif  // INC_4D_MPI_OBSERVABLES_MPI_H
