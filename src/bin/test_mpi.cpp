#include <mpi.h>

#include <print>

#include "../gauge/GaugeField.h"
#include "../geometry/Geometry.h"
#include "../mpi/MpiTopology.h"
#include "../mpi/HalosExchange.h"
#include "../observables/observables.h"

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) std::cout << "Geometry creation..." << std::endl;
    Geometry geo(4, 6);
    if (rank == 0) std::print("MPI topology creation...\n");
    mpi::MpiTopology topo(2);
    if (rank == 0) std::print("Gauge field creation...\n");
    GaugeField field(geo);
    if (rank == 0) std::print("Gauge field ok\n");
    
    std::mt19937_64 rng(123 + rank);
    field.hot_start(geo, rng);
    if (rank == 0) std::print("Exchanging halos...\n");
    mpi::exchange::exchange_halos_cascade(field, geo, topo);
    if (rank == 0) std::print("Exchanging halos ok\n");

    std::cout << "Rank " << rank << " link: " << field.view_link_const(geo.index(1, 1, 1, 3), 3) << "\n";
    if (rank == 0) std::cout << "Plaquette : " << mean_plaquette_weighted(field, geo, 1, 3) << "\n";


    MPI_Finalize();
    return 0;
}
