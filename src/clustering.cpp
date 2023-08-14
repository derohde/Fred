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
        for (curve_number_t i = 0; i < row.size(); ++i) {
            py::print(*row[i], " ");
        }
        py::print();
    }
}

const Curve& Clustering_Result::get(const curve_number_t i) const {
    return centers[i];
}

Cluster_Assignment& Clustering_Result::get_assignment() const {
    return *assignment;
}

void Clustering_Result::set(const curve_number_t i, const Curve &curve) {
    centers[i] = curve;
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
    if (Config::verbosity > 0) py::print("Clustering Result: computing assignment");
    assignment = std::make_unique<Cluster_Assignment>(*this, in, distance_func);
    if (consecutive_call and in.size() == distances.size()) {
        for (curve_number_t i = 0; i < in.size(); ++i) assignment->operator[](_nearest_center(i, in, simplifications, center_indices, distances, distance_func)).push_back(i);
    } else {
        if (Config::use_distance_matrix) distances = Distance_Matrix(in.size(), centers.size());
        
        Curve_Numbers ncenter_indices = Curve_Numbers(centers.size());
        
        for (curve_number_t i = 0; i < centers.size(); ++i) 
            ncenter_indices[i] = i;

        for (curve_number_t i = 0; i < in.size(); ++i) {
            assignment->operator[](_nearest_center(i, in, centers, ncenter_indices, distances, distance_func)).push_back(i);
        }
    }
}

void Clustering_Result::set_center_indices(const Curve_Numbers &pcenter_indices) {
    center_indices = pcenter_indices;
}

py::list Clustering_Result::compute_center_enclosing_balls(const Curves &in, const bool consecutive_call) {
    if (Config::verbosity > 1) py::print("Clustering Result: computing enclosing balls");
    
    py::list result;
    
    compute_assignment(in, consecutive_call);
    
    std::vector<std::vector<Points>> center_matching_points;
    
    for (curve_number_t i = 0; i < size(); ++i) {
        center_matching_points.push_back(std::vector<Points>(get(i).complexity(), get(i).dimensions()));
    }
    
    for (curve_number_t i = 0; i < size(); ++i) {
        if (Config::verbosity > 2) py::print("Clustering Result: computing points for center ", i);
        py::list center_list;
        for (curve_number_t j = 0; j < (*assignment)[i].size(); ++j) {
            Points tpoints(get(i).dimensions());
            switch (distance_func) {
                case 0:
                    if (Config::use_distance_matrix) tpoints = Frechet::Continuous::vertices_matching_points(get(i), in[(*assignment)[i][j]], dynamic_cast<const Frechet::Continuous::Distance&>(*distances[(*assignment)[i][j]][i]));
                    else {
                        Frechet::Continuous::Distance curr_dist = Frechet::Continuous::distance(get(i), in[(*assignment)[i][j]]);
                        tpoints = Frechet::Continuous::vertices_matching_points(get(i), in[(*assignment)[i][j]], curr_dist);
                    }
                    break;
                case 1:
                    break;
                case 2:
                    if (Config::use_distance_matrix) tpoints = Dynamic_Time_Warping::Discrete::vertices_matching_points(get(i), in[(*assignment)[i][j]], dynamic_cast<const Dynamic_Time_Warping::Discrete::Distance&>(*distances[(*assignment)[i][j]][i]));
                    else {
                        Dynamic_Time_Warping::Discrete::Distance curr_dist = Dynamic_Time_Warping::Discrete::distance(get(i), in[(*assignment)[i][j]]);
                        tpoints = Dynamic_Time_Warping::Discrete::vertices_matching_points(get(i), in[(*assignment)[i][j]], curr_dist);
                    }
                    break;
                default:
                    py::print("not implemented!");
            }
            for (curve_size_t k = 0; k < get(i).complexity(); ++k) { 
                center_matching_points[i][k].push_back(tpoints[k]);
            }
        }
        for (curve_size_t k = 0; k < get(i).complexity(); ++k) {
            py::list b_r;
            switch (distance_func) {
                case 0:
                    {
                        const auto bs = bounding_sphere(center_matching_points[i][k]);
                        b_r.append(bs.first.as_ndarray());
                        b_r.append(bs.second);
                        center_list.append(b_r);
                    }
                    break;
                case 1:
                    break;
                case 2:
                    {
                        const auto mean = center_matching_points[i][k].centroid();
                        b_r.append(mean.as_ndarray());
                        b_r.append(0);
                        center_list.append(b_r);
                    }
                    break;
                default:
                    py::print("not implemented!");
            }
        }
        
        
        result.append(center_list);
    }
    
    return result;
}

curve_number_t Cluster_Assignment::count(const curve_number_t i) const {
    return operator[](i).size();
}

curve_number_t Cluster_Assignment::get(const curve_number_t i, const curve_number_t j) const {
    return operator[](i)[j];
}

distance_t Cluster_Assignment::distance(const curve_number_t i, const curve_number_t j) const {
    if (Config::use_distance_matrix) return distances[i][j]->value;
    else {
        switch (distance_func) {
            case 0:
                return Frechet::Continuous::distance(clustering_result.get(i), assignment_curves[(*this)[i][j]]).value;
            case 1:
                return Frechet::Discrete::distance(clustering_result.get(i), assignment_curves[(*this)[i][j]]).value;
            case 2:
                return Dynamic_Time_Warping::Discrete::distance(clustering_result.get(i), assignment_curves[(*this)[i][j]]).value;
            default:
                return std::numeric_limits<distance_t>::signaling_NaN();
        }
    }
}

Clustering_Result kl_cluster(const curve_number_t num_centers, const curve_size_t ell, const Curves &in, unsigned int local_search = 0,
                             const bool median = false, const bool consecutive_call = false, const bool random_start_center = true, const bool fast_simplification = false, const unsigned int distance_func = 0) {
    
    const auto start = std::clock();
    Clustering_Result result(distance_func);
    
    if (in.empty()) return result;
    
    std::size_t memory_distance_matrix = std::pow(in.size(), 2) * (sizeof(distance_t) + 3 * sizeof(double) + ell * in.get_m() * 2 * sizeof(curve_number_t) + sizeof(std::size_t)), 
        memory_available = .666 * Config::available_memory;
        
    if (memory_distance_matrix > memory_available and Config::use_distance_matrix == true) {
        py::print("KL_CLUST: WARNING distance preprocessing requires more memory (", memory_distance_matrix * 1e-9, "GB) than available (", memory_available * 1e-9, "GB), consecutive_call will NOT be available");
        Config::use_distance_matrix = false;
    }

    if (not consecutive_call) {
        if (Config::use_distance_matrix) {
            if (Config::verbosity > 0) py::print("KL_CLUST: allocating ", in.size(), " x ", in.size(), " distance_matrix");
            distances = Distance_Matrix(in.size(), in.size());
        }
        if (Config::verbosity > 0) py::print("KL_CLUST: allocating space for ", in.size(), " simplifications, each of complexity ", ell);
        simplifications = Curves(in.size(), ell, in.dimensions());
    } else if (Config::use_distance_matrix) {
        if (distances.empty()) {
            py::print("WARNING: consecutive_call is used wrongly");
            if (Config::verbosity > 0) py::print("KL_CLUST: allocating ", in.size(), " x ", in.size(), " distance_matrix");
            distances = Distance_Matrix(in.size(), in.size());
            if (Config::verbosity > 0) py::print("KL_CLUST: allocating space for ", in.size(), " simplifications, each of complexity ", ell);
            simplifications = Curves(in.size(), ell, in.dimensions());
        }
        if (distances.size() != in.size()) {
            py::print("WARNING: you have tried to use 'consecutive_call = true' with different input; ignoring!");
            if (Config::verbosity > 0) py::print("KL_CLUST: allocating ", in.size(), " x ", in.size(), " distance_matrix");
            distances = Distance_Matrix(in.size(), in.size());
            if (Config::verbosity > 0) py::print("KL_CLUST: allocating space for ", in.size(), " simplifications, each of complexity ", ell);
            simplifications = Curves(in.size(), ell, in.dimensions());
        }
    }

    Curve_Numbers centers;
    
    const auto simplify = [&](const curve_number_t i) {
        switch (distance_func) {
            case 0:
                if (fast_simplification) {
                    if (Config::verbosity > 0) py::print("KL_CLUST: computing approximate vertex restricted minimum error simplification");
                    auto simplified_curve = Frechet::Continuous::Simplification::approximate_minimum_error_simplification(const_cast<Curve&>(in[i]), ell);
                    simplified_curve.set_name("Simplification of " + in[i].get_name());
                    return simplified_curve;
                } else {
                    if (Config::verbosity > 0) py::print("KL_CLUST: computing exact vertex restricted minimum error simplification");
                    Frechet::Continuous::Simplification::Subcurve_Shortcut_Graph graph(const_cast<Curve&>(in[i]));
                    auto simplified_curve = graph.minimum_error_simplification(ell);
                    simplified_curve.set_name("Simplification of " + in[i].get_name());
                    return simplified_curve;
                }
                break;
            default:
                return in[i];
        }
    };

    if (Config::verbosity > 0) py::print("KL_CLUST: computing first center");
    if (random_start_center) {
        Random::Uniform_Random_Generator<parameter_t> ugen;
        const curve_number_t r =  std::floor(simplifications.size() * ugen.get());
        if (simplifications[r].empty()) {
            if (Config::verbosity > 0) py::print("KL_CLUST: computing simplification of curve ", r);
            simplifications[r] = simplify(r);
        }
        centers.push_back(r);
    } else {
        if (simplifications[0].empty()) {
            if (Config::verbosity > 0) py::print("KL_CLUST: computing simplification of curve 0");
            simplifications[0] = simplify(0);
            centers.push_back(0);
        }
    }
    if (Config::verbosity > 0) py::print("KL_CLUST: first center is ", centers[0]);
    
    distance_t curr_curve_cost, curr_maxdist = 0;
    curve_number_t curr_maxcurve = 0;

    if (Config::verbosity > 0) py::print("KL_CLUST: computing remaining centers");
    {
        // remaining centers
        for (curve_number_t i = 1; i < num_centers; ++i) {
            
            curr_maxdist = 0;
            curr_maxcurve = 0;
            {
            
                if (Config::verbosity > 0) py::print("KL_CLUST: computing new center");
                // all curves
                for (curve_number_t j = 0; j < in.size(); ++j) {
                    
                    curr_curve_cost = _curve_cost(j, in, simplifications, centers, distances, distance_func);
                    
                    if (curr_curve_cost > curr_maxdist) {
                        curr_maxdist = curr_curve_cost;
                        curr_maxcurve = j;
                    }
                    
                }
                if (Config::verbosity > 0) py::print("KL_CLUST: center ", i + 1, " is curve ", curr_maxcurve);
                if (Config::verbosity > 0) py::print("KL_CLUST: current cost is ", curr_maxdist);
                
                if (simplifications[curr_maxcurve].empty()) {
                    if (Config::verbosity > 0) py::print("KL_CLUST: computing simplification of ", curr_maxcurve);
                    simplifications[curr_maxcurve] = simplify(curr_maxcurve);
                }
                centers.push_back(curr_maxcurve);
            }   
        }
    }
    
    if (Config::verbosity > 0) py::print("KL_CLUST: k-center cost is ", curr_maxdist);
    
    if (local_search > 0) {
        
        Curve_Numbers curr_centers = centers;
        distance_t cost = curr_maxdist, curr_cost = cost;
        
        if (Config::verbosity > 0) py::print("KL_CLUST: starting local search for k-center objective for ", local_search, " iterations");
        
        for (unsigned int k = 0; k < local_search; ++k) {
        
            if (Config::verbosity > 0) py::print("KL_CLUST: k-center local search iteration ", k + 1);
        
            for (curve_number_t i = 0; i < num_centers; ++i) {
                
                for (curve_number_t j = 0; j < simplifications.size(); ++j) {
                    
                    if (std::find(curr_centers.begin(), curr_centers.end(), j) != curr_centers.end()) continue;
                    
                    if (Config::verbosity > 0) py::print("KL_CLUST: substituting curve ", curr_centers[i]," for curve ", j," as center");
                    // swap
                    if (simplifications[j].empty()) {
                        if (Config::verbosity > 0) py::print("KL_CLUST: computing simplification of curve ", j);
                        simplifications[j] = simplify(j);
                    }
                    curr_centers[i] = j;
                    // new cost
                    if (Config::verbosity > 0) py::print("KL_CLUST: updating k-center cost");
                    curr_cost = _center_cost_max(in, simplifications, curr_centers, distances, distance_func);
                    // check if improvement is done
                    if (curr_cost < cost) {
                        if (Config::verbosity > 0) py::print("KL_CLUST: cost improves to ", curr_cost);
                        cost = curr_cost;
                        centers = curr_centers;
                    } else {
                        if (Config::verbosity > 0) py::print("KL_CLUST: cost does not improve");
                    }
                }
            }
        }
    }
    
    if (median) {
        
        if (Config::verbosity > 0) py::print("KL_CLUST: computing k-median cost");
        distance_t cost = _center_cost_sum(in, simplifications, centers, distances, distance_func), approxcost = cost, curr_cost = cost;
        if (Config::verbosity > 0) py::print("KL_CLUST: k-median cost is ", cost);
        distance_t gamma = 1/(10 * num_centers);
        bool found = true;
        Curve_Numbers curr_centers = centers;
        
        if (Config::verbosity > 0) py::print("KL_CLUST: starting k-median local search");
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
                                        
                    if (Config::verbosity > 0) py::print("KL_CLUST: substituting curve ", curr_centers[i]," for curve ", j," as center");
                    // swap
                    if (simplifications[j].empty()) {
                        if (Config::verbosity > 0) py::print("KL_CLUST: computing simplification of curve ", j);
                        simplifications[j] = simplify(j);
                    }
                    curr_centers[i] = j;
                    // new cost
                    if (Config::verbosity > 0) py::print("KL_CLUST: updating k-median cost");
                    curr_cost = _center_cost_sum(in, simplifications, curr_centers, distances, distance_func);
                    // check if improvement is done
                    if (curr_cost < cost - gamma * approxcost) {
                        if (Config::verbosity > 0) py::print("KL_CLUST: cost improves to ", curr_cost);
                        cost = curr_cost;
                        centers = curr_centers;
                        found = true;
                    } else {
                        if (Config::verbosity > 0) py::print("KL_CLUST: cost does not improve");
                    }
                }
            }
        }
        curr_maxdist = cost;
    }

    Curves simpl_centers(centers.size(), ell, in.dimensions());
    
    for (curve_number_t i = 0; i < centers.size(); ++i) simpl_centers[i] = simplifications[centers[i]];
    
    const auto end = std::clock();
    result.centers = simpl_centers;
    result.set_center_indices(centers);
    result.value = curr_maxdist;
    result.running_time = (end - start) / CLOCKS_PER_SEC;
    return result;
}

Clustering_Result kl_center(const curve_number_t num_centers, const curve_size_t ell, const Curves &in, unsigned int local_search, const bool consecutive_call, const bool random_start_center, const bool fast_simplification, const unsigned int distance_func) {
    return kl_cluster(num_centers, ell, in, local_search, false, consecutive_call, random_start_center, fast_simplification, distance_func);
}

Clustering_Result kl_median(const curve_number_t num_centers, const curve_size_t ell, const Curves &in, const bool consecutive_call, const bool fast_simplification, const unsigned int distance_func) {
    return kl_cluster(num_centers, ell, in, 0, true, consecutive_call, true, fast_simplification, distance_func);
}

}
