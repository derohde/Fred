/*
Copyright 2020-2021 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <unordered_map>

#include "config.hpp"
#include "types.hpp"
#include "random.hpp"
#include "curve.hpp"
#include "frechet.hpp"
#include "dynamic_time_warping.hpp"
#include "simplification.hpp"

namespace Clustering {
    
struct Distance_Matrix : public std::vector<std::vector<distance_t>> {
    Distance_Matrix() {}
    Distance_Matrix(const curve_number_t n, const curve_number_t m) : std::vector<std::vector<distance_t>>(n, std::vector<distance_t>(m, -1.0)) {}
    
    void print() const;
};

extern Distance_Matrix distances;
extern Curves simplifications;

class Cluster_Assignment : public std::vector<std::vector<curve_number_t>> {
public:
    inline Cluster_Assignment(const curve_number_t k = 0) : std::vector<std::vector<curve_number_t>>(k, std::vector<curve_number_t>()) {}
    
    inline curve_number_t count(const curve_number_t i) const {
        return operator[](i).size();
    }
    
    inline curve_number_t get(const curve_number_t i, const curve_number_t j) const {
        return operator[](i)[j];
    }
    
};

inline distance_t _cheap_dist(const curve_number_t i, const curve_number_t j, const Curves &in, const Curves &simplified_in, Distance_Matrix &distances) {
    if (distances[i][j] < 0) {
        const auto dist = Frechet::Continuous::distance(in[i], simplified_in[j]);
        distances[i][j] = dist.value;
    }
    return distances[i][j];
}

inline curve_number_t _nearest_center(const curve_number_t i, const Curves &in, const Curves &simplified_in, const std::vector<curve_number_t> &centers, Distance_Matrix &distances) {
    const auto infty = std::numeric_limits<distance_t>::infinity();
    // cost for curve is infinity
    auto min_cost = infty;
    curve_number_t nearest = 0;
    
    // except there is a center with smaller cost, then choose the one with smallest cost
    for (curve_number_t j = 0; j < centers.size(); ++j) {
        if (_cheap_dist(i, centers[j], in, simplified_in, distances) < min_cost) {
            min_cost = _cheap_dist(i, centers[j], in, simplified_in, distances);
            nearest = j;
        }
    }
    return nearest;
}

inline distance_t _curve_cost(const curve_number_t i, const Curves &in, const Curves &simplified_in, const std::vector<curve_number_t> &centers, Distance_Matrix &distances) {
    return _cheap_dist(i, centers[_nearest_center(i, in, simplified_in, centers, distances)], in, simplified_in, distances);
}

inline distance_t _center_cost_sum(const Curves &in, const Curves &simplified_in, const std::vector<curve_number_t> &centers, Distance_Matrix &distances) {
    distance_t cost = 0;
    
    // for all curves
    for (curve_number_t i = 0; i < in.size(); ++i) {
        const auto min_cost_elem = _curve_cost(i, in, simplified_in, centers, distances);
        cost += min_cost_elem;
    }
    return cost;
}

struct Clustering_Result {
    Curves centers;
    distance_t value;
    double running_time;
    Cluster_Assignment assignment;
    
    inline Curve& get(const curve_number_t i) {
        return centers[i];
    }
    inline curve_number_t size() const {
        return centers.size();
    }
    
    inline Curves::const_iterator cbegin() const {
        return centers.cbegin();
    }
    
    inline Curves::const_iterator cend() const {
        return centers.cend();
    }
    
    inline void compute_assignment(const Curves &in, const bool consecutive_call = false) {
        assignment = Cluster_Assignment(centers.size());
        if (consecutive_call and in.size() == distances.size()) {
            for (curve_number_t i = 0; i < in.size(); ++i) assignment[_nearest_center(i, in, simplifications, center_indices, distances)].push_back(i);
        } else {
            auto ndistances = Distance_Matrix(in.size(), centers.size());
            auto ncenter_indices = std::vector<curve_number_t>(centers.size());
            std::iota(ncenter_indices.begin(), ncenter_indices.end(), 0);
            for (curve_number_t i = 0; i < in.size(); ++i) assignment[_nearest_center(i, in, centers, ncenter_indices, ndistances)].push_back(i);
        }
    }
    
    inline void set_center_indices(const std::vector<curve_number_t> &pcenter_indices) {
        center_indices = pcenter_indices;
    }

private:
    std::vector<curve_number_t> center_indices;
};

Clustering_Result kl_cluster(const curve_number_t, const curve_size_t, const Curves &, const bool, const bool, const bool, const bool);

Clustering_Result kl_center(const curve_number_t, const curve_size_t, const Curves &, const bool = false, const bool = true, const bool = false);

Clustering_Result kl_median(const curve_number_t, const curve_size_t, const Curves &, const bool = false, const bool = false);

}
