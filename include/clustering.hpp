/*
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <unordered_map>
#include <chrono>

#include "random.hpp"
#include "curve.hpp"
#include "frechet.hpp"
#include "dynamic_time_warping.hpp"
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

inline Cluster_Assignment _cluster_assignment(const Curves &in, const Curves &simplified_in, const std::vector<curve_number_t> &centers, Distance_Matrix &distances) {
    const auto k = centers.size();
    Cluster_Assignment result(centers.size());
    
    if (k == 0) return result;
        
    for (curve_number_t i = 0; i < in.size(); ++i) result[_nearest_center(i, in, simplified_in, centers, distances)].push_back(i);
    
    return result;  
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
    
    inline void compute_assignment(const Curves &in) {
        Distance_Matrix distances(in.size(), centers.size());
        std::vector<curve_number_t> center_indices = std::vector<curve_number_t>(centers.size(), 0);
        for (curve_size_t i = 1; i < centers.size(); ++i) center_indices[i] = i;
        assignment = _cluster_assignment(in, centers, center_indices, distances);
    }
};

Clustering_Result kl_center(const curve_number_t num_centers, const curve_size_t ell, const Curves &in, Distance_Matrix &distances, const bool local_search = false, const Curves &center_domain = Curves(), const bool random_start_center = true) {
    
    const auto start = std::chrono::high_resolution_clock::now();
    Clustering_Result result;
    
    if (in.empty()) return result;
        
    std::vector<curve_number_t> centers;
    Curves &simplified_in = const_cast<Curves&>(center_domain);
    bool self_simplify = false;
    
    if (center_domain.empty()) {
        self_simplify = true;
        Curves simplified_in_self(in.number(), ell, in.dimensions());
        simplified_in = simplified_in_self;
    }
        
    if (random_start_center) {
        Random::Uniform_Random_Generator<double> ugen;
        const curve_number_t r =  std::floor(simplified_in.size() * ugen.get());
        if (self_simplify) {
            Simplification::Subcurve_Shortcut_Graph graph(const_cast<Curve&>(in[r]));
            auto simplified_curve = graph.weak_minimum_error_simplification(ell);
            simplified_curve.set_name("Simplification of " + in[r].get_name());
            simplified_in[r] = simplified_curve;
        }
        centers.push_back(r);
        
    } else {
        if (self_simplify) {
            Simplification::Subcurve_Shortcut_Graph graph(const_cast<Curve&>(in[0]));
            auto simplified_curve = graph.weak_minimum_error_simplification(ell);
            simplified_curve.set_name("Simplification of " + in[0].get_name());
            simplified_in[0] = simplified_curve;
        }
        centers.push_back(0);
    }
    
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
                
                if (self_simplify and simplified_in[curr_maxcurve].empty()) {
                    Simplification::Subcurve_Shortcut_Graph graph(const_cast<Curve&>(in[curr_maxcurve]));
                    auto simplified_curve = graph.weak_minimum_error_simplification(ell);
                    simplified_curve.set_name("Simplification of " + in[curr_maxcurve].get_name());
                    simplified_in[curr_maxcurve] = simplified_curve;
                }
                centers.push_back(curr_maxcurve);
            }   
        }
    }
    
    if (local_search) {
        
        distance_t cost = _center_cost_sum(in, simplified_in, centers, distances);
        distance_t approxcost = cost;
        distance_t curr_cost = cost;
        distance_t gamma = 1/(10 * num_centers);
        bool found = true;
        auto curr_centers = centers;
        
        // try to improve current solution
        while (found) {
            found = false;
            
            // go through all centers
            for (curve_number_t i = 0; i < num_centers; ++i) {
                curr_centers = centers;
                
                // check if there is a better center among all other curves
                for (curve_number_t j = 0; j < simplified_in.size(); ++j) {
                    // continue if curve is already part of center set
                    if (std::find(curr_centers.begin(), curr_centers.end(), j) != curr_centers.end()) continue;
                    
                    // swap
                    if (self_simplify and simplified_in[j].empty()) {
                        Simplification::Subcurve_Shortcut_Graph graph(const_cast<Curve&>(in[j]));
                        auto simplified_curve = graph.weak_minimum_error_simplification(ell);
                        simplified_curve.set_name("Simplification of " + in[j].get_name());
                        simplified_in[j] = simplified_curve;
                    }
                    curr_centers[i] = j;
                    // new cost
                    curr_cost = _center_cost_sum(in, simplified_in, curr_centers, distances);
                    // check if improvement is done
                    if (curr_cost < cost - gamma * approxcost) {
                        cost = curr_cost;
                        centers = curr_centers;
                        found = true;
                    }
                }
            }
        }
        curr_maxdist = cost;
    }

    Curves simpl_centers;
    for (const auto center: centers) simpl_centers.push_back(simplified_in[center]);
    
    auto end = std::chrono::high_resolution_clock::now();
    result.centers = simpl_centers;
    result.value = curr_maxdist;
    result.running_time = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    return result;
}

Clustering_Result kl_median(const curve_number_t num_centers, const curve_size_t ell, const Curves &in, Distance_Matrix &distances, const Curves &center_domain = Curves()) {
    return kl_center(num_centers, ell, in, distances, true, center_domain, false);
}

Clustering_Result one_median_sampling(const curve_size_t ell, const Curves &in, const double epsilon, const Curves &center_domain = Curves()) {
    const auto start = std::chrono::high_resolution_clock::now();
    Clustering_Result result;
    std::vector<curve_number_t> centers;
    Curves &simplified_in = const_cast<Curves&>(center_domain);
    
    if (center_domain.empty()) {
        Curves simplified_in_self(in.number(), ell, in.dimensions());
        
        for (curve_number_t i = 0; i < in.size(); ++i) {
            Simplification::Subcurve_Shortcut_Graph graph(const_cast<Curve&>(in[i]));
            auto simplified_curve = graph.weak_minimum_error_simplification(ell);
            simplified_curve.set_name("Simplification of " + in[i].get_name());
            simplified_in_self[i] = simplified_curve;
        }
        simplified_in = simplified_in_self;
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
    
    auto end = std::chrono::high_resolution_clock::now();
    result.centers.push_back(simplified_in[centers[0]]);
    result.value = _center_cost_sum(in, simplified_in, centers, distances);
    result.running_time = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    return result;
}

Clustering_Result one_median_exhaustive(const curve_size_t ell, const Curves &in, const Curves &center_domain = Curves()) {
    const auto start = std::chrono::high_resolution_clock::now();
    Clustering_Result result;
    std::vector<curve_number_t> centers;
    const Curves &simplified_in = center_domain;
    
    if (center_domain.empty()) {
        Curves simplified_in_self(in.number(), ell, in.dimensions());
        
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
    
    auto end = std::chrono::high_resolution_clock::now();
    result.centers.push_back(simplified_in[centers[0]]);
    result.value = best_objective_value;
    result.running_time = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    return result;
}

Clustering_Result two_two_dtw_one_two_median(const Curves &in, const bool with_assignment = false) {
    const auto start = std::chrono::high_resolution_clock::now();
    Clustering_Result result;
    
    const auto n = in.size();
    
    std::vector<std::vector<bool>> markings = std::vector<std::vector<bool>>(n, std::vector<bool>(in.get_m(), false));
    std::vector<curve_size_t> svert = std::vector<curve_size_t>(n, 0), evert = std::vector<curve_size_t>(n, 0);
    
    Points S1(in.dimensions()), S2(in.dimensions());
    Point mu1(in.dimensions()), mu2(in.dimensions());
    
    for (curve_number_t i = 0; i < n; ++i) {
        S1.push_back(in[i][svert[i]]);
        markings[i][svert[i]] = true;
        ++svert[i];
        evert[i] = in[i].size() - 1;
        S2.push_back(in[i][evert[i]]);
        markings[i][evert[i]] = true;
        --evert[i];
    }
    
    mu1 = S1.centroid();
    mu2 = S2.centroid();
        
    const auto infty = std::numeric_limits<distance_t>::infinity();
    
    bool done = false;
    distance_t d1 = infty, d2 = infty, dist = infty;
    curve_number_t c1 = 0, c2 = 0;
    
    while (not done) {
        d1 = infty;
        d2 = infty;
        
        for (curve_size_t i = 0; i < in.size(); ++i) {
            if (not markings[i][svert[i]]) {
                dist = in[i][svert[i]].dist_sqr(mu1);
                if (dist < d1) {
                    d1 = dist;
                    c1 = i;
                }
            }
            if (not markings[i][evert[i]]) {
                dist = in[i][evert[i]].dist_sqr(mu2);
                if (dist < d2) {
                    d2 = dist;
                    c2 = i;
                }
            }
        }
        
        if (d1 < d2) {
            //std::cout << "S1 add " << c1 << "." << svert[c1] << std::endl;
            S1.push_back(in[c1][svert[c1]]);
            markings[c1][svert[c1]] = true;
            ++svert[c1];
            mu1 = S1.centroid();
            done = false;
        }
        else if (d2 < infty) {
            //std::cout << "S2 add " << c2 << "." << evert[c2] << std::endl;
            S2.push_back(in[c2][evert[c2]]);
            markings[c2][evert[c2]] = true;
            --evert[c2];
            mu2 = S2.centroid();
            done = false;
        } else done = true;
        
    }
    
    Curve center_curve(mu1.dimensions(), "center curve");
    center_curve.push_back(mu1);
    center_curve.push_back(mu2);
    
    distance_t cost = 0;
    
    for (const auto &p : S1) cost += p.dist(mu1);
    for (const auto &p : S2) cost += p.dist(mu2);
    
    auto end = std::chrono::high_resolution_clock::now();
    result.centers.push_back(center_curve);
    result.value = cost;
    result.running_time = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    return result;
}

Clustering_Result two_two_dtw_one_two_median_exact(const Curves &in, const bool with_assignment = false) {
    const auto start = std::chrono::high_resolution_clock::now();
    Clustering_Result result;
    Curve best_center(in.dimensions());
    const auto infty = std::numeric_limits<distance_t>::infinity();
    
    curve_size_t n = 1;
    std::vector<curve_size_t> pointers = std::vector<curve_size_t>(in.size(), 0), 
                                divisors = std::vector<curve_size_t>(in.size(), 0);
    distance_t best = infty, cost = 0;
    Points S1(in.dimensions()), S2(in.dimensions());
    
    for (curve_number_t i = 0; i < in.size(); ++i) {
        n *= in[i].complexity() - 1;
        if (i == 0) divisors[i] = in[0].complexity() - 1;
        else if (i == 1) divisors[i] = in[0].complexity() - 1;
        else divisors[i] = divisors[i-1] * (in[i].complexity() - 1);
    }
    
    const auto onepercent = n / 100;
    
    int currperc = 0;
    
    for (curve_size_t i = 0; i < n; ++i) {
        
        if (onepercent > 0) {
            if (i / onepercent > currperc) {
                currperc = i / onepercent;
                std::cout << currperc << "% done" << std::endl;
            }
        }
        
        pointers[0] = i % divisors[0];
        for (curve_number_t j = 1; j < in.size(); ++j) {
            pointers[j] = (i / divisors[j]) % (in[j].complexity() - 1);
        }
        
        S1.clear();
        S2.clear();
        
        for (curve_number_t j = 0; j < in.size(); ++j) {
            for (curve_size_t k = 0; k < in[j].complexity(); ++ k) {
                if (k <= pointers[j]) {
                    S1.push_back(in[j][k]);
                    if (k == in[j].complexity() - 1) std::cerr << "error!!" << std::endl;
                }
                else S2.push_back(in[j][k]);
            }
        }
        
        auto mu1 = S1.centroid();
        auto mu2 = S2.centroid();
        
        Curve center_curve(mu1.dimensions(), "optimal center curve");
        center_curve.push_back(mu1);
        center_curve.push_back(mu2);
        
        cost = 0;
        
        for (curve_number_t j = 0; j < in.size(); ++j) {
            const auto dist = Dynamic_Time_Warping::Discrete::distance(center_curve, in[j]);
            cost += dist.value;
        }
        
        if (cost < best) {
            best = cost;
            best_center = center_curve;
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    result.centers.push_back(best_center);
    result.value = best;
    result.running_time = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    return result;
}


}
