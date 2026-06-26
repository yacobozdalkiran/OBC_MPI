#ifndef OBC_MPI_OVERRELAXATION_H
#define OBC_MPI_OVERRELAXATION_H

#include "../gauge/GaugeField.h"
#include "../su3/utils.h"

namespace overrelaxation {
    // Calculates the overrelaxation step for a single SU(2) subgroup block w
    SU2q su2_step(const SU2q &w);

    // Performs an overrelaxation update on a single link
    void hit(GaugeField &field, const Geometry &geo, size_t site, int mu, SU3 &A);

    // Full overrelaxation sweep (sequential/single-thread)
    void sweep(GaugeField &field, const Geometry &geo);
}

namespace mpi::overrelaxationcb {
    // Performs an overrelaxation hit on a single link (used by checkerboard)
    void hit(GaugeField &field, const Geometry &geo, size_t site, int mu, SU3 &A);

    // Performs an overrelaxation sweep on links of the specified parity
    void sweep(GaugeField &field, const Geometry &geo, site_parity update_parity);

    // Runs a specified number of overrelaxation sweeps using checkerboard updates
    void sample(GaugeField &field, const Geometry &geo, int N_sweeps);
}

#endif // OBC_MPI_OVERRELAXATION_H
