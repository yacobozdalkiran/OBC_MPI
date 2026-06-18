#include "Geometry.h"

#include <cstdint>

// Initializes the precomputed tables neighbors, frozen and list_staples
Geometry::Geometry(int T_, int L_) {
    T = T_;
    L_int = L_;
    V_int = L_int * L_int * L_int * T;
    L_ext = L_int + 2;
    V_ext = L_ext * L_ext * L_ext * T;
    neighbors.resize((V_ext) * 8, SIZE_MAX);
    for (int t = 0; t < T; t++) {
        for (int z = 0; z < L_ext; z++) {
            for (int y = 0; y < L_ext; y++) {
                for (int x = 0; x < L_ext; x++) {
                    size_t site = index(x, y, z, t);
                    // Direction X (mu=0)
                    if (x + 1 < L_ext) neighbors[index_neigh(site, 0, up)] = index(x + 1, y, z, t);
                    if (x - 1 >= 0) neighbors[index_neigh(site, 0, down)] = index(x - 1, y, z, t);
                    // Direction Y (mu=1)
                    if (y + 1 < L_ext) neighbors[index_neigh(site, 1, up)] = index(x, y + 1, z, t);
                    if (y - 1 >= 0) neighbors[index_neigh(site, 1, down)] = index(x, y - 1, z, t);
                    // Direction Z (mu=2)
                    if (z + 1 < L_ext) neighbors[index_neigh(site, 2, up)] = index(x, y, z + 1, t);
                    if (z - 1 >= 0) neighbors[index_neigh(site, 2, down)] = index(x, y, z - 1, t);
                    // Direction T (mu=3)
                    if (t + 1 < T) neighbors[index_neigh(site, 3, up)] = index(x, y, z, t + 1);
                    if (t - 1 >= 0) neighbors[index_neigh(site, 3, down)] = index(x, y, z, t - 1);
                }
            }
        }
    }

    frozen.resize(V_ext * 4, false);
    // Frozen links are those that step out of the lattice core or belong to halos
    for (int t = 0; t < T; t++) {
        for (int z = 0; z < L_ext; z++) {
            for (int y = 0; y < L_ext; y++) {
                for (int x = 0; x < L_ext; x++) {
                    bool link_is_frozen = false;
                    // Halos are always frozen
                    if (x == 0 or y == 0 or z == 0) link_is_frozen = true;
                    if (x == L_ext - 1 or y == L_ext - 1 or z == L_ext - 1)
                        link_is_frozen = true;
                    for (int mu = 0; mu < 4; mu++) {
                        // Parallelization
                        if (x == 1 and mu != 0) link_is_frozen = true;
                        if (y == 1 and mu != 1) link_is_frozen = true;
                        if (z == 1 and mu != 2) link_is_frozen = true;
                        // if (t == 1 and mu != 3) link_is_frozen = true; (because no
                        // parallelization in t)
                        size_t i = index(x, y, z, t);
                        frozen[index_frozen(i, mu)] = link_is_frozen;
                    }
                }
            }
        }
    }

    parity.resize(V_ext);
    for (int t = 0; t < T; t++) {
        for (int z = 0; z < L_ext; z++) {
            for (int y = 0; y < L_ext; y++) {
                for (int x = 0; x < L_ext; x++) {
                    size_t i = index(x, y, z, t);
                    parity[i] = ((x + y + z + t) % 2 == 0) ? EV : OD;
                }
            }
        }
    }

    links_staples.resize(V_ext * 4 * 6 * 3, std::make_pair(SIZE_MAX, -1));
    // 4 links per site, 6 staples per link, 3 links per staple
    for (int t = 0; t < T; t++) {
        for (int z = 1; z < L_ext - 1; z++) {
            for (int y = 1; y < L_ext - 1; y++) {
                for (int x = 1; x < L_ext - 1; x++) {
                    size_t site = index(x, y, z, t);  // x
                    for (int mu = 0; mu < 4; mu++) {
                        int j = 0;
                        for (int nu = 0; nu < 4; nu++) {
                            if (nu == mu) continue;
                            if (mu == 3 and t == T - 1) continue;

                            // Forward staple
                            size_t xmu = get_neigh(site, mu, up);  // x+mu
                            size_t xnu = get_neigh(site, nu, up);  // x+nu
                            // Backward staple
                            size_t xmunu = get_neigh(xmu, nu, down);  // x+mu-nu
                            size_t xmnu = get_neigh(site, nu, down);  // x-nu
                            // We fill only the valid staples (ie not containing a link fixed to Id)
                            if (not(t == T - 1 and mu == 3) and not(t == T - 1 and nu == 3)) {
                                links_staples[index_staples(site, mu, j, 0)] = {xmu, nu};
                                links_staples[index_staples(site, mu, j, 1)] = {xnu, mu};
                                links_staples[index_staples(site, mu, j, 2)] = {site, nu};
                            }

                            if (not(t == 0 and nu == 3) and not(t == T - 1 and mu == 3)) {
                                links_staples[index_staples(site, mu, j + 1, 0)] = {xmunu, nu};
                                links_staples[index_staples(site, mu, j + 1, 1)] = {xmnu, mu};
                                links_staples[index_staples(site, mu, j + 1, 2)] = {xmnu, nu};
                            }
                            j += 2;
                        }
                    }
                }
            }
        }
    }

    staple_valid.resize(V_ext * 4 * 6, false);
    // 4 links per site, 6 staples per link
    for (int t = 0; t < T; t++) {
        for (int z = 1; z < L_ext - 1; z++) {
            for (int y = 1; y < L_ext - 1; y++) {
                for (int x = 1; x < L_ext - 1; x++) {
                    size_t site = index(x, y, z, t);  // x
                    for (int mu = 0; mu < 4; mu++) {
                        int j = 0;
                        for (int nu = 0; nu < 4; nu++) {
                            if (nu == mu) continue;
                            if (mu == 3 and t == T - 1) continue;
                            // We mark the valid staples for OBC (ie not containing a link fixed to
                            // Id)
                            if (not(t == T - 1 and mu == 3) and not(t == T - 1 and nu == 3)) {
                                staple_valid[index_staple_valid(site, mu, j)] = true;
                            }

                            if (not(t == 0 and nu == 3) and not(t == T - 1 and mu == 3)) {
                                staple_valid[index_staple_valid(site, mu, j + 1)] = true;
                            }
                            j += 2;
                        }
                    }
                }
            }
        }
    }

    staple_coeff.resize(V_ext * 4 * 6, 0);
    // 4 links per site, 6 staples per link
    for (int t = 0; t < T; t++) {
        for (int z = 1; z < L_ext - 1; z++) {
            for (int y = 1; y < L_ext - 1; y++) {
                for (int x = 1; x < L_ext - 1; x++) {
                    size_t site = index(x, y, z, t);
                    for (int mu = 0; mu < 4; mu++) {
                        int j = 0;
                        for (int nu = 0; nu < 4; nu++) {
                            if (nu == mu) continue;
                            if (mu == 3 and t == T - 1) continue;
                            // We mark the valid staples (ie not containing a link fixed to Id)
                            // Forward staple
                            if (not(t == T - 1 and mu == 3) and not(t == T - 1 and nu == 3)) {
                                // Spatial plaquette of a border spatial link
                                if ((t == T - 1 and mu != 3 and nu != 3) or
                                    (t == 0 and mu != 3 and nu != 3))
                                    staple_coeff[index_staple_valid(site, mu, j)] = 0.5;
                                else
                                    staple_coeff[index_staple_valid(site, mu, j)] = 1.0;
                            }
                            // Backward staple
                            if (not(t == 0 and nu == 3) and not(t == T - 1 and mu == 3)) {
                                // Spatial plaquette of a border spatial link
                                if ((t == T - 1 and mu != 3 and nu != 3) or
                                    (t == 0 and mu != 3 and nu != 3))
                                    staple_coeff[index_staple_valid(site, mu, j + 1)] = 0.5;
                                else
                                    staple_coeff[index_staple_valid(site, mu, j + 1)] = 1.0;
                            }
                            j += 2;
                        }
                    }
                }
            }
        }
    }
    // --- Optimization: Flattening valid staples ---
    fwd_start.assign(V_ext * 4 + 1, 0);
    bwd_start.assign(V_ext * 4 + 1, 0);

    for (size_t site = 0; site < V_ext; ++site) {
        for (int mu = 0; mu < 4; ++mu) {
            size_t link_idx = site * 4 + mu;
            fwd_start[link_idx + 1] = fwd_start[link_idx];
            bwd_start[link_idx + 1] = bwd_start[link_idx];

            for (int j = 0; j < 6; j += 2) {
                // Forward staple
                if (is_staple_valid(site, mu, j)) {
                    auto l0 = get_link_staple(site, mu, j, 0);
                    auto l1 = get_link_staple(site, mu, j, 1);
                    auto l2 = get_link_staple(site, mu, j, 2);
                    fwd_staples_opt.push_back(
                        {(l0.first * 4 + l0.second) * 9, (l1.first * 4 + l1.second) * 9,
                         (l2.first * 4 + l2.second) * 9, get_staple_coeff(site, mu, j), j});
                    fwd_start[link_idx + 1]++;
                }

                // Backward staple
                if (is_staple_valid(site, mu, j + 1)) {
                    auto l3 = get_link_staple(site, mu, j + 1, 0);
                    auto l4 = get_link_staple(site, mu, j + 1, 1);
                    auto l5 = get_link_staple(site, mu, j + 1, 2);
                    bwd_staples_opt.push_back(
                        {(l3.first * 4 + l3.second) * 9, (l4.first * 4 + l4.second) * 9,
                         (l5.first * 4 + l5.second) * 9, get_staple_coeff(site, mu, j + 1), j + 1});
                    bwd_start[link_idx + 1]++;
                }
            }
        }
    }
}
