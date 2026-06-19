#include <mpi.h>
#include <print>
#include <iostream>
#include <cmath>

#include "../gauge/GaugeField.h"
#include "../geometry/Geometry.h"
#include "../mpi/MpiTopology.h"
#include "../mpi/HalosExchange.h"
#include "../io/ildg.h"

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    mpi::MpiTopology topo(2);
    Geometry geo(8, 4); // T = 8, L = 4
    GaugeField field(geo);

    // Populate the field with deterministic values (including temporal slice patterns)
    for (int t = 0; t < geo.T; t++) {
        for (int z = 1; z <= geo.L_int; z++) {
            for (int y = 1; y <= geo.L_int; y++) {
                for (int x = 1; x <= geo.L_int; x++) {
                    size_t site = geo.index(x, y, z, t);
                    for (int mu = 0; mu < 4; mu++) {
                        SU3 U = SU3::Zero();
                        // Put recognizable, coordinate-dependent values
                        double val = 1.0 + t * 1000.0 + z * 100.0 + y * 10.0 + x + mu * 0.1;
                        U(0, 0) = Complex(val, -val * 0.5);
                        U(1, 1) = Complex(val * 2.0, val * 0.25);
                        U(2, 2) = Complex(val * 0.5, -val * 2.0);
                        field.view_link(site, mu) = U;
                    }
                }
            }
        }
    }

    std::string filename = "test_config";
    std::string dirpath = "./test_io";

    if (topo.rank == 0) {
        std::print("Saving configuration to ILDG format...\n");
    }
    save_ildg_clime(filename, dirpath, field, geo, topo);

    if (topo.rank == 0) {
        std::print("Loading configuration back from ILDG format...\n");
    }
    GaugeField field_read(geo);
    read_ildg_clime(filename, dirpath, field_read, geo, topo);

    // Verify
    bool success = true;
    for (int t = 0; t < geo.T; t++) {
        for (int z = 1; z <= geo.L_int; z++) {
            for (int y = 1; y <= geo.L_int; y++) {
                for (int x = 1; x <= geo.L_int; x++) {
                    size_t site = geo.index(x, y, z, t);
                    for (int mu = 0; mu < 4; mu++) {
                        auto original = field.view_link_const(site, mu);
                        auto loaded = field_read.view_link_const(site, mu);
                        for (int i = 0; i < 3; i++) {
                            for (int j = 0; j < 3; j++) {
                                double diff_re = std::abs(original(i, j).real() - loaded(i, j).real());
                                double diff_im = std::abs(original(i, j).imag() - loaded(i, j).imag());
                                if (diff_re > 1e-12 || diff_im > 1e-12) {
                                    std::print("Rank {} mismatch at t={}, z={}, y={}, x={}, mu={}, link({}, {}):\n",
                                               topo.rank, t, z, y, x, mu, i, j);
                                    std::print("  Original: ({}, {})\n", original(i, j).real(), original(i, j).imag());
                                    std::print("  Loaded:   ({}, {})\n", loaded(i, j).real(), loaded(i, j).imag());
                                    success = false;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    int local_ok = success ? 1 : 0;
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (topo.rank == 0) {
        if (global_ok == topo.size) {
            std::print("ILDG SAVE/READ VERIFICATION SUCCESSFUL!\n");
        } else {
            std::print("ILDG SAVE/READ VERIFICATION FAILED!\n");
        }
    }

    MPI_Finalize();
    return global_ok == topo.size ? 0 : 1;
}
