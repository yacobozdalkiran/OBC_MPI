//
// Created by ozdalkiran-l on 1/14/26.
//

#ifndef ECMC_MPI_PARAMS_H
#define ECMC_MPI_PARAMS_H

#include <string>
struct ECMCParams {
    double beta = 6.0;
    double param_theta_sample = 100;
    double param_theta_refresh = 50;
    bool poisson = false;
};

struct HbParams {
    double beta = 6.0;
    int N_hits = 1;
    int N_sweeps = 1;
};

// Test purpose
struct RunParamsECMC {
    int T = 8;
    int L = 4;
    bool cold_start = true;
    int seed = 123;
    ECMCParams ecmc_params{};  // Params of the ECMC for each even/odd update
    int N_samples = 10;
    bool topo = true;
    int N_plaquette = 2;  // Measure plaquette every N_shift_plaquette_shift
    int T_min = 0;
    int T_max = T - 1;
    std::string run_name = "c";
    std::string run_dir = "data";
    int save_each = 2;  // save confs/measures/seed each
};

// Test purpose
struct RunParamsHb {
    int T = 8;
    int L = 6;
    bool cold_start = true;
    int seed = 123;
    HbParams hp{};  // Hb params for each even/odd update
    int N_samples = 10;
    bool topo = true;
    int N_plaquette = 2;  // Measure plaquette every N_shift_plaquette_shift
    int T_min = 0;
    int T_max = T - 1;
    std::string run_name = "c";
    std::string run_dir = "data";
    int save_each = 2;  // save confs/measures/seed each
};

// Main
struct RunParamsECB {
    int T = 8;
    int L_core = 6;
    int n_core_dims = 2;
    bool cold_start = true;

    int seed = 123;
    int N_shift = 2;  // Number of shifts

    int N_ov_sweep = 2;  // Nb ov sweeps
    int N_ov_hit = 1;    // Nb ov hits

    ECMCParams ecmc_params{};   // Params of the ECMC for each even/odd update
                                //
    int N_shift_plaquette = 2;  // Measure plaquette every N_shift_plaquette_shift
    int T_min = 0;
    int T_max = T - 1;

    bool topo = true;
    int N_shift_topo = 10;  // Measure topo charge every N_shift_topo shift
    int N_steps_gf = 10;
    int N_rk_steps = 40;

    std::string run_name = "c";
    std::string run_dir = "data";
    int save_each_shifts = 2;  // save confs/measures/seed each
};

// Main
struct RunParamsHbCB {
    int T = 8;
    int L_core = 6;
    int n_core_dims = 2;
    bool cold_start = true;

    int seed = 123;
    int N_shift = 2;  // Number of shifts

    int N_ov_sweep = 2;  // Nb ov sweeps
    int N_ov_hit = 1;    // Nb ov hits

    HbParams hp{};  // Hb params for each even/odd update

    int N_shift_plaquette = 2;  // Measure plaquette every N_shift_plaquette_shift
    int T_min = 0;
    int T_max = T - 1;

    bool topo = true;
    int N_shift_topo = 10;  // Measure topo charge every N_shift_topo shift
    int N_steps_gf = 10;
    int N_rk_steps = 40;

    std::string run_name = "c";
    std::string run_dir = "data";
    int save_each_shifts = 2;  // save confs/measures/seed each
};

#endif  // ECMC_MPI_PARAMS_H
