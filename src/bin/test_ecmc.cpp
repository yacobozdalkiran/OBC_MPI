#include <print>

#include "../ecmc/ecmc.h"
#include "../gauge/GaugeField.h"
#include "../io/params.h"
#include "../mpi/HalosExchange.h"
#include "../mpi/Shift.h"
#include "../observables/observables_mpi.h"

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    // Objects initialization

    // MPI
    mpi::MpiTopology topo(2);
    // Rng
    int seed = topo.rank * 1000;
    std::mt19937_64 rng(seed);
    // Geometry
    Geometry geo(10, 6);
    // Field
    GaugeField field(geo);
    // Halos for shifts
    HalosShift h(geo);

    // ECMC
    ECMCParams ep{.beta = 6.0,
                  .param_theta_sample = 6000,
                  .param_theta_refresh = 800,
                  .poisson = false,
                  .epsilon_set = 0.15};
    LocalChainState chain{};
    chain.initialized = false;
    Distributions d{ep};

    // Field hot start
    field.cold_start(geo);
    mpi::exchange::exchange_halos_cascade(field, geo, topo);

    double plaquette{};
    plaquette = mpi::observables::mean_plaquette_global(field, geo, topo, 3, 6);
    if (topo.rank == 0) std::print("Initial plaquette : {}\n", plaquette);

    int N_shifts = 100;
    for (int shifts = 0; shifts < N_shifts; shifts++) {
        mpi::shift::random_shift(field, geo, h, topo, rng);
        mpi::exchange::exchange_halos_cascade(field, geo, topo);
        ecmc::sample_persistant_norev(chain, d, field, geo, ep, rng);
        mpi::exchange::exchange_halos_cascade(field, geo, topo);
        plaquette = mpi::observables::mean_plaquette_global(field, geo, topo, 0, 9);
        if (topo.rank == 0) std::print("Shift {}, plaquette : {}\n", shifts, plaquette);
    }

    MPI_Finalize();
}
