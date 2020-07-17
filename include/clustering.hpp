/*
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <unordered_map>

#include <boost/chrono/include.hpp>

#include "random.hpp"
#include "curve.hpp"
#include "frechet.hpp"
#include "simplification.hpp"

namespace Clustering {
    
    
struct Distance_Matrix : public std::vector<std::vector<distance_t>> {
    Distance_Matrix() {}
    Distance_Matrix(const curve_number_t n, const curve_number_t m) : std::vector<std::vector<distance_t>>(n, std::vector<distance_t>(m, -1.0)) {}
    
    void print() const {
        for (const auto &row : *this) {
            for (const auto elem : row) {
                std::cout << elem << " ";
            }
            std::cout << std::endl;
        }
    }
};


class Cluster_Assignment : public std::unordered_map<curve_number_t, std::vector<curve_number_t>> {
public:
    
    inline const std::vector<curve_number_t>& operator[](const curve_number_t i) const {
        return std::unordered_map<curve_number_t, std::vector<curve_number_t>>::at(i);
    }
    
    inline std::vector<curve_number_t>& operator[](const curve_number_t i) {
        return std::unordered_map<curve_number_t, std::vector<curve_number_t>>::at(i);
    }
    
    inline curve_number_t count(const curve_number_t i) const {
        return operator[](i).size();
    }
    
    inline curve_number_t get(const curve_number_t i, const curve_number_t j) const {
        return operator[](i)[j];
    }
    
};

struct Clustering_Result {
    Curves centers;
    distance_t value;
    double running_time;
    Cluster_Assignment assignment;
    
    inline Curve get(const curve_number_t i) const {
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
};


inline void _cheap_dist(const curve_number_t i, const curve_number_t j, const Curves &in, const Curves &simplified_in, Distance_Matrix &distances) {
    if (distances[i][j] < 0) {
        const auto dist = Frechet::Continuous::distance(in[i], simplified_in[j]);
        distances[i][j] = dist.value;
    }
}

inline curve_number_t _nearest_center(const curve_number_t i, const Curves &in, const Curves &simplified_in, const std::vector<curve_number_t> &centers, Distance_Matrix &distances) {
    const auto infty = std::numeric_limits<distance_t>::infinity();
    // cost for curve is infinity
    auto min_cost_elem = infty;
    curve_number_t nearest = 0;
    
    // except there is a center with smaller cost, then choose the one with smallest cost
    for (curve_number_t j = 0; j < centers.size(); ++j) {
        _cheap_dist(i, centers[j], in, simplified_in, distances);
        if (distances[i][centers[j]] < min_cost_elem) {
            min_cost_elem = distances[i][centers[j]];
            nearest = j;
        }
    }
    return nearest;
}

inline distance_t _curve_cost(const curve_number_t i, const Curves &in, const Curves &simplified_in, const std::vector<curve_number_t> &centers, Distance_Matrix &distances) {
    const auto infty = std::numeric_limits<distance_t>::infinity();
    // cost for curve is infinity
    auto min_cost_elem = infty;
    curve_number_t nearest = 0;
    
    // except there is a center with smaller cost, then choose the one with smallest cost
    for (curve_number_t j = 0; j < centers.size(); ++j) {
        _cheap_dist(i, centers[j], in, simplified_in, distances);
        if (distances[i][centers[j]] < min_cost_elem) {
            min_cost_elem = distances[i][centers[j]];
            nearest = j;
        }
    }
    return min_cost_elem;
}

inline distance_t _center_cost_sum(const Curves &in, const Curves &simplified_in, const std::vector<curve_number_t> &centers, Distance_Matrix &distances) {
    distance_t cost = 0.0;
    
    // for all curves
    for (curve_number_t i = 0; i < in.size(); ++i) {
        const auto min_cost_elem = _curve_cost(i, in, simplified_in, centers, distances);
        cost += min_cost_elem;
    }
    return cost;
}

inline Cluster_Assignment _cluster_assignment(const Curves &in, const Curves &simplified_in, const std::vector<curve_number_t> &centers, Distance_Matrix &distances) {
    Cluster_Assignment result;
    const auto k = centers.size();
    
    if (k == 0) return result;
    
    for (curve_number_t i = 0; i < k; ++i) result.emplace(i, std::vector<curve_number_t>());
    
    for (curve_number_t i = 0; i < in.size(); ++i) result[_nearest_center(i, in, simplified_in, centers, distances)].push_back(i);
    
    return result;  
}

Clustering_Result gonzalez(const curve_number_t num_centers, const curve_size_t ell, const Curves &in, Distance_Matrix &distances, const bool arya = false, const bool with_assignment = false, 
                           const Curves &center_domain = Curves(), const bool random_start_center = true) {
    
    const auto start = boost::chrono::process_real_cpu_clock::now();
    Clustering_Result result;
    
    if (in.empty()) return result;
        
    std::vector<curve_number_t> centers;
    const Curves &simplified_in = center_domain;
    
    if (center_domain.empty()) {
        Curves simplified_in_self(in.number(), ell);
        
        for (curve_number_t i = 0; i < in.size(); ++i) {
            Simplification::Subcurve_Shortcut_Graph graph(const_cast<Curve&>(in[i]));
            auto simplified_curve = graph.weak_minimum_error_simplification(ell);
            simplified_curve.set_name("Simplification of " + in[i].get_name());
            simplified_in_self[i] = simplified_curve;
        }
        const_cast<Curves&>(simplified_in) = simplified_in_self;
    }
    
    const auto n = in.size();
    
    if (random_start_center) {
        
        Random::Uniform_Random_Generator<double> ugen;
        const curve_number_t r =  std::floor(n * ugen.get());
        centers.push_back(r);
        
    } else centers.push_back(0);
    
    distance_t curr_maxdist = 0;
    curve_number_t curr_maxcurve = 0;
    distance_t curr_curve_cost;
    
    if (distances.empty()) distances = Distance_Matrix(in.size(), simplified_in.size());
    
    {
        // remaining centers
        for (curve_number_t i = 1; i < num_centers; ++i) {
            
            curr_maxdist = 0;
            curr_maxcurve = 0;
            {
            
                // all curves
                for (curve_number_t j = 0; j < in.size(); ++j) {
                    
                    curr_curve_cost = _curve_cost(j, in, simplified_in, centers, distances);
                    
                    if (curr_curve_cost > curr_maxdist) {
                        curr_maxdist = curr_curve_cost;
                        curr_maxcurve = j;
                    }
                    
                }
                #if DEBUG
                std::cout << "found center no. " << i+1 << std::endl;
                #endif
                
                centers.push_back(curr_maxcurve);
            }   
        }
    }
    
    if (arya) {
        
        auto cost = _center_cost_sum(in, simplified_in, centers, distances);
        auto approxcost = cost;
        auto gamma = 1/(3 * num_centers * in.size());
        auto found = false;
        
        // try to improve current solution
        while (true) {
            found = false;
            
            // go through all centers
            for (curve_number_t i = 0; i < num_centers; ++i) {
                auto curr_centers = centers;
                
                // check if there is a better center among all other curves
                for (curve_number_t j = 0; j < in.size(); ++j) {
                    // continue if curve is already part of center set
                    if (std::find(curr_centers.begin(), curr_centers.end(), j) != curr_centers.end()) continue;
                    
                    // swap
                    curr_centers[i] = j;
                    // new cost
                    auto curr_cost = _center_cost_sum(in, simplified_in, curr_centers, distances);
                    // check if improvement is done
                    if (cost - gamma * approxcost > curr_cost) {
                        cost = curr_cost;
                        centers = curr_centers;
                        found = true;
                    }
                }
            }
            if (not found) break;
        }
        curr_maxdist = cost;
    }
    
    if (with_assignment) {
        result.assignment = _cluster_assignment(in, simplified_in, centers, distances);
    }
    
    Curves simpl_centers;
    for (const auto center: centers) simpl_centers.push_back(simplified_in[center]);
    
    auto end = boost::chrono::process_real_cpu_clock::now();
    result.centers = simpl_centers;
    result.value = curr_maxdist;
    result.running_time = (end-start).count() / 1000000000.0;
    return result;
}

Clustering_Result arya(const curve_number_t num_centers, const curve_size_t ell, const Curves &in, Distance_Matrix &distances, 
                       const bool with_assignment = false, const Curves &center_domain = Curves(), const bool random_start_center = true) {
    return gonzalez(num_centers, ell, in, distances, true, with_assignment, center_domain, random_start_center);
}

Clustering_Result one_median_sampling(const curve_size_t ell, const Curves &in, const double epsilon, const bool with_assignment = false, const Curves &center_domain = Curves()) {
    const auto start = boost::chrono::process_real_cpu_clock::now();
    Clustering_Result result;
    std::vector<curve_number_t> centers;
    const Curves &simplified_in = center_domain;
    
    if (center_domain.empty()) {
        Curves simplified_in_self(in.number(), ell);
        
        for (curve_number_t i = 0; i < in.size(); ++i) {
            Simplification::Subcurve_Shortcut_Graph graph(const_cast<Curve&>(in[i]));
            auto simplified_curve = graph.weak_minimum_error_simplification(ell);
            simplified_curve.set_name("Simplification of " + in[i].get_name());
            simplified_in_self[i] = simplified_curve;
        }
        const_cast<Curves&>(simplified_in) = simplified_in_self;
    }
    
    const auto n = in.size();
    
    const auto s = std::ceil(60);
    const auto t = std::ceil(std::log(60)/(epsilon*epsilon));
    
    Random::Uniform_Random_Generator<double> ugen;
    
    const auto candidates = ugen.get(s);
    const auto witnesses = ugen.get(t);
    
    Distance_Matrix distances = Distance_Matrix(in.size(), in.size());
    
    curve_number_t best_candidate = 0;
    distance_t best_objective_value = std::numeric_limits<distance_t>::infinity();
    
    for (curve_number_t i = 0; i < candidates.size(); ++i) {
        
        const curve_number_t candidate = std::floor(candidates[i] * n);
        distance_t objective = 0;
        
        for (curve_number_t j = 0; j < witnesses.size(); ++j) {
            const curve_number_t witness = std::floor(witnesses[j] * n);
            
            _cheap_dist(witness, candidate, in, simplified_in, distances);
            objective += distances[witness][candidate];
        }
        
        if (objective < best_objective_value) {
            best_candidate = candidate;
            best_objective_value = objective;
        }
    }
    centers.push_back(best_candidate);
    
    if (with_assignment) {
        result.assignment = _cluster_assignment(in, simplified_in, centers, distances);
    }
    
    auto end = boost::chrono::process_real_cpu_clock::now();
    result.centers.push_back(simplified_in[centers[0]]);
    result.value = _center_cost_sum(in, simplified_in, centers, distances);
    result.running_time = (end-start).count() / 1000000000.0;
    return result;
}

Clustering_Result one_median_exhaustive(const curve_size_t ell, const Curves &in, const bool with_assignment = false, const Curves &center_domain = Curves()) {
    const auto start = boost::chrono::process_real_cpu_clock::now();
    Clustering_Result result;
    std::vector<curve_number_t> centers;
    const Curves &simplified_in = center_domain;
    
    if (center_domain.empty()) {
        Curves simplified_in_self(in.number(), ell);
        
        for (curve_number_t i = 0; i < in.size(); ++i) {
            Simplification::Subcurve_Shortcut_Graph graph(const_cast<Curve&>(in[i]));
            auto simplified_curve = graph.weak_minimum_error_simplification(ell);
            simplified_curve.set_name("Simplification of " + in[i].get_name());
            simplified_in_self[i] = simplified_curve;
        }
        const_cast<Curves&>(simplified_in) = simplified_in_self;
    }
    
    const auto n = in.size();
        
    Distance_Matrix distances = Distance_Matrix(in.size(), in.size());
    
    curve_number_t best_candidate = 0;
    distance_t best_objective_value = std::numeric_limits<distance_t>::infinity();
    
    for (curve_number_t j = 0; j < in.size(); ++j) {
        
        distance_t objective = 0;
        
        for (curve_number_t i = 0; i < in.size(); ++i) {
            _cheap_dist(i, j, in, simplified_in, distances);
            objective += distances[i][j];
        }
        
        if (objective < best_objective_value) {
            best_candidate = j;
            best_objective_value = objective;
        }
    }
    centers.push_back(best_candidate);
    
    if (with_assignment) {
        result.assignment = _cluster_assignment(in, simplified_in, centers, distances);
    }
    
    auto end = boost::chrono::process_real_cpu_clock::now();
    result.centers.push_back(simplified_in[centers[0]]);
    result.value = best_objective_value;
    result.running_time = (end-start).count() / 1000000000.0;
    return result;
}

}
