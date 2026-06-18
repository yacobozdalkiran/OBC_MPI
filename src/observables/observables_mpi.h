#ifndef INC_4D_MPI_OBSERVABLES_MPI_H
#define INC_4D_MPI_OBSERVABLES_MPI_H

#include "../gauge/GaugeField.h"
#include "../mpi/MpiTopology.h"
#include "observables.h"

namespace mpi::observables {
double mean_plaquette_global(GaugeField& field, const Geometry& geo, MpiTopology& topo, int T_min,
                             int T_max);
SU3 clover_site(const GaugeField& field, const Geometry& geo, size_t site, int mu, int nu);
int levi_civita(int mu, int nu, int rho, int sigma);
}  // namespace mpi::observables

#endif  // INC_4D_MPI_OBSERVABLES_MPI_H
