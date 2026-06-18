//
// Created by ozdalkiran-l on 1/28/26.
//

#ifndef ECMC_MPI_HEATBATH_MPI_H
#define ECMC_MPI_HEATBATH_MPI_H

#include "../gauge/GaugeField.h"
#include "../io/params.h"
#include "../su3/utils.h"

namespace mpi::heatbathcb {
void hit(GaugeField& field, const Geometry& geo, size_t site, int mu, double beta, SU3& A,
         std::mt19937_64& rng);
void sweep(GaugeField& field, const Geometry& geo, double beta, int N_hits,
           std::vector<std::mt19937_64>& rng, site_parity update_parity);
std::vector<double> samples(GaugeField& field, const Geometry& geo, const HbParams& params,
                            std::vector<std::mt19937_64>& rng);
}  // namespace mpi::heatbathcb

#endif  // ECMC_MPI_HEATBATH_MPI_H
