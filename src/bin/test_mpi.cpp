#include <mpi.h>

#include <print>

#include "../gauge/GaugeField.h"
#include "../geometry/Geometry.h"
#include "../mpi/MpiTopology.h"
#include "../mpi/HalosExchange.h"
#include "../observables/observables_mpi.h"

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
    
    // --- Verification of halo exchange with hardcoded values ---
    if (rank == 0) std::print("Populating gauge field with rank/coordinate patterns...\n");
    for (int t = 0; t < geo.T; t++) {
        for (int z = 1; z <= geo.L_int; z++) {
            for (int y = 1; y <= geo.L_int; y++) {
                for (int x = 1; x <= geo.L_int; x++) {
                    size_t site = geo.index(x, y, z, t);
                    for (int mu = 0; mu < 4; mu++) {
                        SU3 U = SU3::Zero();
                        U(0, 0) = Complex(rank * 1000 + x * 100 + y * 10 + z, 0.0);
                        field.view_link(site, mu) = U;
                    }
                }
            }
        }
    }

    if (rank == 0) std::print("Exchanging halos...\n");
    mpi::exchange::exchange_halos_cascade(field, geo, topo);
    if (rank == 0) std::print("Exchanging halos ok. Validating...\n");

    bool success = true;
    
    // Check 1: minus-X halo
    {
        size_t site = geo.index(0, 1, 1, 0);
        double val = field.view_link_const(site, 0)(0, 0).real();
        double expected = topo.x0 * 1000 + geo.L_int * 100 + 1 * 10 + 1;
        if (std::abs(val - expected) > 1e-5) {
            std::print("Rank {} Error: x=0 halo got {}, expected {}\n", rank, val, expected);
            success = false;
        }
    }
    
    // Check 2: minus-Y halo
    {
        size_t site = geo.index(1, 0, 1, 0);
        double val = field.view_link_const(site, 0)(0, 0).real();
        double expected = topo.y0 * 1000 + 1 * 100 + geo.L_int * 10 + 1;
        if (std::abs(val - expected) > 1e-5) {
            std::print("Rank {} Error: y=0 halo got {}, expected {}\n", rank, val, expected);
            success = false;
        }
    }

    // Check 3: minus-Z halo
    {
        size_t site = geo.index(1, 1, 0, 0);
        double val = field.view_link_const(site, 0)(0, 0).real();
        double expected = topo.z0 * 1000 + 1 * 100 + 1 * 10 + geo.L_int;
        if (std::abs(val - expected) > 1e-5) {
            std::print("Rank {} Error: z=0 halo got {}, expected {}\n", rank, val, expected);
            success = false;
        }
    }

    // Check 4: diagonal X-Y halo (cascade validation)
    {
        int xy_coords[3] = {
            (topo.local_coords[0] - 1 + 2) % 2,
            (topo.local_coords[1] - 1 + 2) % 2,
            topo.local_coords[2]
        };
        int xy_rank;
        MPI_Cart_rank(topo.cart_comm, xy_coords, &xy_rank);
        size_t site = geo.index(0, 0, 1, 0);
        double val = field.view_link_const(site, 0)(0, 0).real();
        double expected = xy_rank * 1000 + geo.L_int * 100 + geo.L_int * 10 + 1;
        if (std::abs(val - expected) > 1e-5) {
            std::print("Rank {} Error: x=0,y=0 diagonal halo got {}, expected {}\n", rank, val, expected);
            success = false;
        }
    }

    int local_ok = success ? 1 : 0;
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    if (rank == 0) {
        if (global_ok == size) {
            std::print("ALL HALO EXCHANGE TESTS PASSED SUCCESSFULLY!\n\n");
        } else {
            std::print("SOME HALO EXCHANGE TESTS FAILED! (passed count: {}/{})\n\n", global_ok, size);
        }
    }

    // Re-initialize with hot start and exchange halos
    std::mt19937_64 rng(123 + rank);
    field.hot_start(geo, rng);
    field.cold_start(geo);
    mpi::exchange::exchange_halos_cascade(field, geo, topo);

    double p = mpi::observables::mean_plaquette_global(field, geo, topo, 1, 3);
    if (topo.rank == 0) std::print("Plaquette MPI : {}\n", p);
    MPI_Finalize();
    return 0;
}
