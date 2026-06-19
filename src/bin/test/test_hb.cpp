#include <print>

#include "../../heatbath/heatbath_mpi.h"
#include "../../gauge/GaugeField.h"
#include "../../io/params.h"
#include "../../mpi/HalosExchange.h"
#include "../../mpi/Shift.h"
#include "../../observables/observables_mpi.h"

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    // Objects initialization

    // MPI
    mpi::MpiTopology topo(2);
    // Vector rng for OpenMP
    int n_threads = omp_get_max_threads();
    std::vector<std::mt19937_64> rng(n_threads);
    for (int i = 0; i < n_threads; i++) {
        // On multiplie le rank par 1000 pour éviter tout recouvrement de séquence
        rng[i].seed(123 + topo.rank * 1000 + i);
    }
    // Geometry
    Geometry geo(8, 4);
    // Field
    GaugeField field(geo);
    // Halos for shifts
    HalosShift h(geo);

    // Heatbath
    HbParams hp{
        .beta=6.0,
        .N_hits=1,
        .N_sweeps=1,
    };

    // Field hot start
    field.cold_start(geo);
    mpi::exchange::exchange_halos_cascade(field, geo, topo);

    double plaquette{};
    plaquette = mpi::observables::mean_plaquette_global(field, geo, topo, 3, 6);
    if (topo.rank == 0) std::print("Initial plaquette : {}\n", plaquette);

    int N_shifts = 100;
    for (int shifts = 0; shifts < N_shifts; shifts++) {
        mpi::shift::random_shift(field, geo, h, topo, rng[0]);
        mpi::exchange::exchange_halos_cascade(field, geo, topo);
        mpi::heatbathcb::sample(field, geo, hp, rng);
        mpi::exchange::exchange_halos_cascade(field, geo, topo);
        plaquette = mpi::observables::mean_plaquette_global(field, geo, topo, 0, 7);
        if (topo.rank == 0) std::print("Shift {}, plaquette : {}\n", shifts, plaquette);
    }

    MPI_Finalize();
}
