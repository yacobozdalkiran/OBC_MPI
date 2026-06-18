#ifndef GEOMETRYCB_H
#define GEOMETRYCB_H

#include <vector>

#include "types.h"

class Geometry {
public:
    int T;
    int L_int;
    int L_ext;
    size_t V_int;
    size_t V_ext;

private:
    std::vector<size_t> neighbors;
    std::vector<bool> frozen;
    std::vector<std::pair<size_t, int>> links_staples;
    std::vector<site_parity> parity;
    std::vector<bool> staple_valid;
    std::vector<double> staple_coeff;

public:
    explicit Geometry(int T_, int L_);

    // Index function for links
    [[nodiscard]] size_t index(int x, int y, int z, int t) const {
        return ((static_cast<size_t>(t) * L_ext + z) * L_ext + y) * L_ext + x;
    }

    // Index of a site in neighbor vector
    [[nodiscard]] static size_t index_neigh(size_t site, int mu, dir d) {
        return site * 8 + mu * 2 + d;
    }

    // Index of a link in is_frozen
    [[nodiscard]] static size_t index_frozen(size_t site, int mu) { 
        return site * 4 + mu; 
    }

    // Index of a link in links_staples
    [[nodiscard]] static size_t index_staples(size_t site, int mu, int i_staple, int i_link) {
        return site * 3 * 6 * 4 + mu * 3 * 6 + i_staple * 3 + i_link;
    }

    // Index of a staple in staple_valid
    [[nodiscard]] static size_t index_staple_valid(size_t site, int mu, int i_staple) {
        return site * 24 + mu * 6 + i_staple;
    }

    // Get a neighbor
    [[nodiscard]] size_t get_neigh(size_t site, int mu, dir d) const {
        return neighbors[index_neigh(site, mu, d)];
    }

    // Returns true if link is frozen
    [[nodiscard]] bool is_frozen(size_t site, int mu) const {
        return frozen[index_frozen(site, mu)];
    }

    // Returns the coords of the link of index i_link of the staple of index i_staple of <site,mu>
    [[nodiscard]] std::pair<size_t, int> get_link_staple(size_t site, int mu, int i_staple,
                                                         int i_link) const {
        return links_staples[index_staples(site, mu, i_staple, i_link)];
    }

    [[nodiscard]] site_parity get_parity(size_t site) const{
        return parity[site];
    }

    //Returns true if the staple i_staple of link site,mu is valid in OBC
    [[nodiscard]] bool is_staple_valid(size_t site, int mu, int i_staple) const {
        return staple_valid[index_staple_valid(site, mu, i_staple)];
    }

    // Return the normalization coefficient of the staple i_staple of link site,mu and 0 if invalid staple
    [[nodiscard]] double get_staple_coeff(size_t site, int mu, int i_staple) const {
        return staple_coeff[index_staple_valid(site, mu, i_staple)];
    }

    // Optimized structure to store pre-calculated memory offsets for staple links
    // off0, off1, off2 are the raw indices into the GaugeField::links vector
    struct OptimizedStaple {
        size_t off0, off1, off2;
        double coeff;
        int j_local; // Index 0-5 identifying the staple slot
    };

    // Flattened storage for only valid staples to avoid branches during simulation
    std::vector<OptimizedStaple> fwd_staples_opt;
    std::vector<OptimizedStaple> bwd_staples_opt;

    // Start index in the *_staples_opt vectors for each link (site * 4 + mu)
    // The range of staples for a link is [start[i], start[i+1])
    std::vector<size_t> fwd_start;
    std::vector<size_t> bwd_start;
};

#endif
