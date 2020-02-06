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

namespace Clustering {

class Centers : public std::vector<curve_size_t> {
public:
    inline const auto get(const curve_size_t i) const {
        return const_cast<Centers*>(this)->operator[](i);
    }
};

class Cluster_Assignment : public std::unordered_map<curve_size_t, std::vector<curve_size_t>> {
public:
    inline const auto count(const curve_size_t i) const {
        return const_cast<Cluster_Assignment*>(this)->operator[](i).size();
    }
    
    inline const auto get(const curve_size_t i, const curve_size_t j) const {
        return const_cast<Cluster_Assignment*>(this)->operator[](i)[j];
    }
    
};

struct Clustering_Result {
    Centers centers;
    distance_t value;
    double running_time;
    Cluster_Assignment assignment;
    
    inline auto get(const curve_size_t i) const {
        return centers.get(i);
    }
    inline auto size() const {
        return centers.size();
    }
};


inline void cheap_dist(const curve_size_t i, const curve_size_t j, const Curves &in, std::vector<std::vector<distance_t>> &distances, 
                        const distance_t eps, const bool round) {
    if (distances[i][j] < 0) {
        const auto dist = Frechet::Continuous::distance(in[i], in[j], eps, round);
        distances[j][i] = dist.value;
        distances[i][j] = dist.value;
    }
}

inline curve_size_t nearest_center(const curve_size_t i, const Curves &in, const Centers &centers, 
        std::vector<std::vector<distance_t>> &distances, const distance_t eps, const bool round) {
    const auto infty = std::numeric_limits<distance_t>::infinity();
    // cost for curve is infinity
    auto min_cost_elem = infty;
    curve_size_t nearest = 0;
    
    // except there is a center with smaller cost, then choose the one with smallest cost
    for (curve_size_t j = 0; j < centers.size(); ++j) {
        cheap_dist(i, centers[j], in, distances, eps, round);
        if (distances[i][centers[j]] < min_cost_elem) {
            min_cost_elem = distances[i][centers[j]];
            nearest = j;
        }
    }
    return nearest;
}

inline auto curve_cost(const curve_size_t i, const Curves &in, const Centers &centers, 
                        std::vector<std::vector<distance_t>> &distances, const distance_t eps, const bool round) {
    const auto infty = std::numeric_limits<distance_t>::infinity();
    // cost for curve is infinity
    auto min_cost_elem = infty;
    
    // except there is a center with smaller cost, then choose the one with smallest cost
    for (curve_size_t j = 0; j < centers.size(); ++j) {
        cheap_dist(i, centers[j], in, distances, eps, round);
        if (distances[i][centers[j]] < min_cost_elem) {
            min_cost_elem = distances[i][centers[j]];
        }
    }
    
    return min_cost_elem;
}

inline auto center_cost_sum(const Curves &in, const Centers &centers, 
                                std::vector<std::vector<distance_t>> &distances, const distance_t eps, const bool round) {
    distance_t cost = 0.0;
    
    // for all curves
    for (curve_size_t i = 0; i < in.size(); ++i) {
        const auto min_cost_elem = curve_cost(i, in, centers, distances, eps, round);
        cost += min_cost_elem;
    }
    return cost;
}

inline Cluster_Assignment cluster_assignment(const Curves &in, const Centers &centers, 
                                                std::vector<std::vector<distance_t>> &distances, const distance_t eps, const bool round) {
    Cluster_Assignment result;
    const auto k = centers.size();
    
    if (k == 0) return result;
    
    for (curve_size_t i = 0; i < k; ++i) result.emplace(i, std::vector<curve_size_t>());
    
    for (curve_size_t i = 0; i < in.size(); ++i) result[nearest_center(i, in, centers, distances, eps, round)].push_back(i);
    
    return result;  
}

Clustering_Result gonzalez(const curve_size_t num_centers, const Curves &in, const distance_t eps, 
                                const bool round = true, const bool arya = false, const bool with_assignment = false) {
    const auto start = boost::chrono::process_real_cpu_clock::now();
    Clustering_Result result;
    
    if (in.empty()) return result;
        
    Centers centers;
    
    const auto n = in.size();
    
    centers.push_back(0);
    
    distance_t curr_maxdist = 0;
    curve_size_t curr_maxcurve = 0;
    
    std::vector<std::vector<distance_t>> distances(in.size(), std::vector<distance_t>(in.size(), -1.0));
    
    // no cost for distances from curves to themselves
    for (curve_size_t i = 0; i < in.size(); ++i) distances[i][i] = 0;
    
    {
        // remaining centers
        for (curve_size_t i = 1; i < num_centers; ++i) {
            
            curr_maxdist = 0;
            curr_maxcurve = 0;
            {
            
                // all curves
                for (curve_size_t j = 0; j < in.size(); ++j) {
                    
                    auto curr_curve_cost = curve_cost(j, in, centers, distances, eps, round);
                    
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
        
        auto cost = center_cost_sum(in, centers, distances, eps, round);
        auto approxcost = cost;
        auto gamma = 1/(3 * num_centers * in.size());
        auto found = false;
        
        // try to improve current solution
        while (true) {
            found = false;
            
            // go through all centers
            for (curve_size_t i = 0; i < num_centers; ++i) {
                auto curr_centers = centers;
                
                // check if there is a better center among all other curves
                for (curve_size_t j = 0; j < in.size(); ++j) {
                    // continue if curve is already part of center set
                    if (std::find(curr_centers.begin(), curr_centers.end(), j) != curr_centers.end()) continue;
                    
                    // swap
                    curr_centers[i] = j;
                    // new cost
                    auto curr_cost = center_cost_sum(in, curr_centers, distances, eps, round);
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
        result.assignment = cluster_assignment(in, centers, distances, eps, round);
    }
    
    auto end = boost::chrono::process_real_cpu_clock::now();
    result.centers = centers;
    result.value = curr_maxdist;
    result.running_time = (end-start).count() / 1000000000.0;
    return result;
}

Clustering_Result arya(const curve_size_t num_centers, const Curves &in, const distance_t eps, const bool round = true, const bool with_assignment = false) {
    return gonzalez(num_centers, in, eps, round, true, with_assignment);
}

Clustering_Result one_median_sampling(const double epsilon, const Curves &in, const distance_t eps, 
                                        const bool round = true, const bool with_assignment = false) {
    const auto start = boost::chrono::process_real_cpu_clock::now();
    Clustering_Result result;
    Centers centers;
    
    const auto n = in.size();
    
    const auto s = std::ceil(60);
    const auto t = std::ceil(std::log(60)/(epsilon*epsilon));
    
    Random::Uniform_Random_Generator<double> ugen;
    
    const auto candidates = ugen.get(s);
    const auto witnesses = ugen.get(t);
    
    std::vector<std::vector<distance_t>> distances(in.size(), 
        std::vector<distance_t>(in.size(), -1.0));
    
    curve_size_t best_candidate = 0;
    distance_t best_objective_value = std::numeric_limits<distance_t>::infinity();
    
    for (curve_size_t i = 0; i < candidates.size(); ++i) {
        
        const curve_size_t candidate = std::floor(candidates[i] * n);
        distance_t objective = 0;
        
        for (curve_size_t j = 0; j < witnesses.size(); ++j) {
            const curve_size_t witness = std::floor(witnesses[j] * n);
            
            cheap_dist(candidate, witness, in, distances, eps, round);
            objective += distances[candidate][witness];
        }
        
        if (objective < best_objective_value) {
            best_candidate = candidate;
            best_objective_value = objective;
        }
    }
    centers.push_back(best_candidate);
    
    if (with_assignment) {
        result.assignment = cluster_assignment(in, centers, distances, eps, round);
    }
    
    auto end = boost::chrono::process_real_cpu_clock::now();
    result.centers = centers;
    result.value = center_cost_sum(in, centers, distances, eps, round);
    result.running_time = (end-start).count() / 1000000000.0;
    return result;
}

Clustering_Result one_median_exhaustive(const Curves &in, const distance_t eps, const bool round = true, const bool with_assignment = false) {
    const auto start = boost::chrono::process_real_cpu_clock::now();
    Clustering_Result result;
    Centers centers;
    
    const auto n = in.size();
        
    std::vector<std::vector<distance_t>> distances(in.size(), 
        std::vector<distance_t>(in.size(), -1.0));
    
    curve_size_t best_candidate = 0;
    distance_t best_objective_value = std::numeric_limits<distance_t>::infinity();
    
    for (curve_size_t i = 0; i < in.size(); ++i) {
        
        distance_t objective = 0;
        
        for (curve_size_t j = 0; j < in.size(); ++j) {
            cheap_dist(i, j, in, distances, eps, round);
            objective += distances[i][j];
        }
        
        if (objective < best_objective_value) {
            best_candidate = i;
            best_objective_value = objective;
        }
    }
    centers.push_back(best_candidate);
    
    if (with_assignment) {
        result.assignment = cluster_assignment(in, centers, distances, eps, round);
    }
    
    auto end = boost::chrono::process_real_cpu_clock::now();
    result.centers = centers;
    result.value = best_objective_value;
    result.running_time = (end-start).count() / 1000000000.0;
    return result;
}

}
