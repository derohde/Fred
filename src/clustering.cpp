/*
Copyright 2021 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include<ctime>

#include "clustering.hpp"

namespace Clustering {

Distance_Matrix distances;
Curves simplifications;

void Distance_Matrix::print() const {
    for (const auto &row : *this) {
        for (const auto elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
}

Curve& Clustering_Result::get(const curve_number_t i) {
    return centers[i];
}

curve_number_t Clustering_Result::size() const {
    return centers.size();
}

Curves::const_iterator Clustering_Result::cbegin() const {
    return centers.cbegin();
}

Curves::const_iterator Clustering_Result::cend() const {
    return centers.cend();
}

void Clustering_Result::compute_assignment(const Curves &in, const bool consecutive_call) {
    assignment = Cluster_Assignment(centers.size());
    if (consecutive_call and in.size() == distances.size()) {
        for (curve_number_t i = 0; i < in.size(); ++i) assignment[_nearest_center(i, in, simplifications, center_indices, distances)].push_back(i);
    } else {
        auto ndistances = Distance_Matrix(in.size(), centers.size());
        auto ncenter_indices = Curve_Numbers(centers.size());
        std::iota(ncenter_indices.begin(), ncenter_indices.end(), 0);
        for (curve_number_t i = 0; i < in.size(); ++i) assignment[_nearest_center(i, in, centers, ncenter_indices, ndistances)].push_back(i);
    }
}

void Clustering_Result::set_center_indices(const Curve_Numbers &pcenter_indices) {
    center_indices = pcenter_indices;
}

curve_number_t Cluster_Assignment::count(const curve_number_t i) const {
    return operator[](i).size();
}

curve_number_t Cluster_Assignment::get(const curve_number_t i, const curve_number_t j) const {
    return operator[](i)[j];
}

Clustering_Result kl_cluster(const curve_number_t num_centers, const curve_size_t ell, const Curves &in, 
                             const bool local_search = false, const bool consecutive_call = false, const bool random_start_center = true, const bool fast_simplification = false) {
    
    const auto start = std::clock();
    Clustering_Result result;
    
    if (in.empty()) return result;

    if (not consecutive_call) {
        distances = Distance_Matrix(in.size(), in.size());
        simplifications = Curves(in.number(), ell, in.dimensions());
    } else {
        if (distances.empty()) {
            distances = Distance_Matrix(in.size(), in.size());
            simplifications = Curves(in.number(), ell, in.dimensions());
        }
        if (distances.size() != in.size()) {
            std::cerr << "WARNING: you have tried to use 'consecutive_call = true' with different input; ignoring!" << std::endl;
            distances = Distance_Matrix(in.size(), in.size());
            simplifications = Curves(in.number(), ell, in.dimensions());
        }
    }

    Curve_Numbers centers;
    
    const auto simplify = [&](const curve_number_t i) {
        if (fast_simplification) {
            if (Config::verbosity > 0) std::cout << "KL_CLUST: computing approximate vertex restricted minimum error simplification" << std::endl;
            auto simplified_curve = Simplification::approximate_minimum_error_simplification(const_cast<Curve&>(in[i]), ell);
            simplified_curve.set_name("Simplification of " + in[i].get_name());
            return simplified_curve;
        } else {
            if (Config::verbosity > 0) std::cout << "KL_CLUST: computing exact vertex restricted minimum error simplification" << std::endl;
            Simplification::Subcurve_Shortcut_Graph graph(const_cast<Curve&>(in[i]));
            auto simplified_curve = graph.minimum_error_simplification(ell);
            simplified_curve.set_name("Simplification of " + in[i].get_name());
            return simplified_curve;
        }
    };

    if (Config::verbosity > 0) std::cout << "KL_CLUST: computing first center" << std::endl;
    if (random_start_center) {
        Random::Uniform_Random_Generator<double> ugen;
        const curve_number_t r =  std::floor(simplifications.size() * ugen.get());
        if (simplifications[r].empty()) {
            if (Config::verbosity > 0) std::cout << "KL_CLUST: computing simplification of curve " << r << std::endl;
            simplifications[r] = simplify(r);
        }
        centers.push_back(r);
    } else {
        if (simplifications[0].empty()) {
            if (Config::verbosity > 0) std::cout << "KL_CLUST: computing simplification of curve 0" << std::endl;
            simplifications[0] = simplify(0);
        }
        centers.push_back(0);
    }
    if (Config::verbosity > 0) std::cout << "KL_CLUST: first center is " << centers[0] << std::endl;
    
    distance_t curr_maxdist = 0;
    curve_number_t curr_maxcurve = 0;
    distance_t curr_curve_cost;

    if (Config::verbosity > 0) std::cout << "KL_CLUST: computing remaining centers" << std::endl;
    {
        // remaining centers
        for (curve_number_t i = 2; i <= num_centers; ++i) {
            
            curr_maxdist = 0;
            curr_maxcurve = 0;
            {
            
                if (Config::verbosity > 0) std::cout << "KL_CLUST: computing new center " << std::endl;
                // all curves
                for (curve_number_t j = 0; j < in.size(); ++j) {
                    
                    curr_curve_cost = _curve_cost(j, in, simplifications, centers, distances);
                    
                    if (curr_curve_cost > curr_maxdist) {
                        curr_maxdist = curr_curve_cost;
                        curr_maxcurve = j;
                    }
                    
                }
                if (Config::verbosity > 0) std::cout << "KL_CLUST: center " << i << " is curve " << curr_maxcurve << std::endl;
                
                if (simplifications[curr_maxcurve].empty()) {
                    if (Config::verbosity > 0) std::cout << "KL_CLUST: computing simplification of " << curr_maxcurve << std::endl;
                    simplifications[curr_maxcurve] = simplify(curr_maxcurve);
                }
                centers.push_back(curr_maxcurve);
            }   
        }
    }
    
    if (local_search) {
        
        if (Config::verbosity > 0) std::cout << "KL_CLUST: computing k-median cost" << std::endl;
        distance_t cost = _center_cost_sum(in, simplifications, centers, distances);
        distance_t approxcost = cost;
        distance_t curr_cost = cost;
        distance_t gamma = 1/(10 * num_centers);
        bool found = true;
        auto curr_centers = centers;
        
        if (Config::verbosity > 0) std::cout << "KL_CLUST: starting local search" << std::endl;
        // try to improve current solution
        while (found) {
            found = false;
            
            // go through all centers
            for (curve_number_t i = 0; i < num_centers; ++i) {
                curr_centers = centers;
                
                // check if there is a better center among all other curves
                for (curve_number_t j = 0; j < simplifications.size(); ++j) {
                    // continue if curve is already part of center set
                    if (std::find(curr_centers.begin(), curr_centers.end(), j) != curr_centers.end()) continue;
                                        
                    if (Config::verbosity > 0) std::cout << "KL_CLUST: substituting curve " << curr_centers[i] << " for curve " << j << " as center" << std::endl;
                    // swap
                    if (simplifications[j].empty()) {
                        if (Config::verbosity > 0) std::cout << "KL_CLUST: computing simplification of curve " << j << std::endl;
                        simplifications[j] = simplify(j);
                    }
                    curr_centers[i] = j;
                    // new cost
                    if (Config::verbosity > 0) std::cout << "KL_CLUST: updating k-median cost" << std::endl;
                    curr_cost = _center_cost_sum(in, simplifications, curr_centers, distances);
                    // check if improvement is done
                    if (curr_cost < cost - gamma * approxcost) {
                        if (Config::verbosity > 0) std::cout << "KL_CLUST: cost did improve" << std::endl;
                        cost = curr_cost;
                        centers = curr_centers;
                        found = true;
                    } else {
                        if (Config::verbosity > 0) std::cout << "KL_CLUST: cost did not improve" << std::endl;
                    }
                }
            }
        }
        curr_maxdist = cost;
    }

    Curves simpl_centers;
    for (const auto center: centers) simpl_centers.push_back(simplifications[center]);
    
    const auto end = std::clock();
    result.centers = simpl_centers;
    result.set_center_indices(centers);
    result.value = curr_maxdist;
    result.running_time = (end - start) / CLOCKS_PER_SEC;
    return result;
}

Clustering_Result kl_center(const curve_number_t num_centers, const curve_size_t ell, const Curves &in, const bool consecutive_call, const bool random_start_center, const bool fast_simplification) {
    return kl_cluster(num_centers, ell, in, false, consecutive_call, random_start_center, fast_simplification);
}

Clustering_Result kl_median(const curve_number_t num_centers, const curve_size_t ell, const Curves &in, const bool consecutive_call, const bool fast_simplification) {
    return kl_cluster(num_centers, ell, in, true, consecutive_call, true, fast_simplification);
}

}
