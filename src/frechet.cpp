/*
Copyright 2020-2021 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <vector>
#include <limits>
#include <ctime>

#include "frechet.hpp"

namespace Frechet {

namespace Continuous {
    
distance_t error = 1;
bool round = true;
    
std::string Distance::repr() const {
    std::stringstream ss;
    ss << value;
    return ss.str();
}

Points vertices_matching_points(const Curve &curve1, const Curve &curve2, const distance_t dist) {
    if (Config::verbosity > 1) py::print("CFD: computing matching points for distance ", dist, " from curve1 of complexity ", curve1.complexity(), " to curve2 of complexity ", curve2.complexity());
    if (Config::verbosity > 2) py::print("CFD: distance between curve1 and curve2 is ", distance(curve1, curve2).value);
    if ((curve1.complexity() < 2) or (curve2.complexity() < 2)) {
        py::print("WARNING: curves must be of at least two points");
        Points result(curve1.dimensions());
        return result;
    }
    
    const auto dist_sqr = dist * dist;
    const curve_size_t n1 = curve1.complexity();
    const curve_size_t n2 = curve2.complexity();

    std::vector<Intervals> free_intervals(n1, Intervals(n2, Interval()));
        
    #pragma omp parallel for collapse(2) if (n1 * n2 > 1000)
    for (curve_size_t i = 1; i < n1; ++i) {
        for (curve_size_t j = 0; j < n2 - 1; ++j) {
            free_intervals[i][j] = curve1[i].ball_intersection_interval(dist_sqr, curve2[j], curve2[j+1]);
        }
    }
    
    if (Config::verbosity > 1) py::print("CFD: free space computed, computing matching");
    
    Points result(n1, curve1.dimensions());
    parameter_t p;
    curve_size_t jj(0);
        
    for (curve_size_t i = 1; i < n1 - 1; ++i) {
        if (Config::verbosity > 1) py::print("CFD: computing matching points for vertex ", i);
        for (curve_size_t j = jj; j < n2 - 1; ++j, p = 0) {
            if (not free_intervals[i][j].empty()) {
                if (j == jj) {
                    p = std::max(p, free_intervals[i][j].begin());
                    break;
                }
                p = free_intervals[i][j].begin();
                jj = j;
                break;
            }
        }
        result[i] = curve2[jj].line_segment_point(curve2[jj+1], p);
        if (Config::verbosity > 1) py::print("CFD: matching vertex ", i, "to ", p, "on segment ", jj, " to ", jj+1, " with distance ", result[i].dist(curve1[i]));
    }
    result[0] = curve2[0];
    result[n1-1] = curve2[n2-1];
    
    return result;
}

Distance distance(const Curve &curve1, const Curve &curve2) {
    if ((curve1.complexity() < 2) or (curve2.complexity() < 2)) {
        py::print("WARNING: comparison possible only for curves of at least two points");
        Distance result;
        result.value = std::numeric_limits<distance_t>::signaling_NaN();
        return result;
    }
    if (curve1.dimensions() != curve2.dimensions()) {
        py::print("WARNING: comparison possible only for curves of equal number of dimensions");
        Distance result;
        result.value = std::numeric_limits<distance_t>::signaling_NaN();
        return result;
    }
    
    const auto start = std::clock();
    if (Config::verbosity > 2) py::print("CFD: computing lower bound");
    const distance_t lb = _projective_lower_bound(curve1, curve2);
    if (Config::verbosity > 2) py::print("CFD: computing upper bound");
    const distance_t ub = _greedy_upper_bound(curve1, curve2);
    const auto end = std::clock();
    
    auto dist = _distance(curve1, curve2, ub, lb);
    dist.time_bounds = (end - start) / CLOCKS_PER_SEC;

    return dist;
}

Distance _distance(const Curve &curve1, const Curve &curve2, distance_t ub, distance_t lb) {
    Distance result;
    const auto start = std::clock();
    
    distance_t split = (ub + lb)/2;
    const distance_t p_error = lb * error / 100 > std::numeric_limits<distance_t>::epsilon() ? lb * error / 100 : std::numeric_limits<distance_t>::epsilon();
    std::size_t number_searches = 0;
    
    if (ub - lb > p_error) {
        if (Config::verbosity > 2) py::print("CFD: binary search using FSD, error = ", p_error);
        
        const auto infty = std::numeric_limits<parameter_t>::infinity();
        std::vector<Parameters> reachable1(curve1.complexity() - 1, Parameters(curve2.complexity(), infty));
        std::vector<Parameters> reachable2(curve1.complexity(), Parameters(curve2.complexity() - 1, infty));
        
        std::vector<Intervals> free_intervals1(curve2.complexity(), Intervals(curve1.complexity(), Interval()));
        std::vector<Intervals> free_intervals2(curve1.complexity(), Intervals(curve2.complexity(), Interval()));

        if (std::isnan(lb) or std::isnan(ub)) {
            result.value = std::numeric_limits<distance_t>::signaling_NaN();
            return result;
        }

        bool isLessThan;
        
        //Binary search over the feasible distances
        while (ub - lb > p_error) {
            ++number_searches;
            split = (ub + lb)/distance_t(2);
            if (split == lb or split == ub) break;
            isLessThan = _less_than_or_equal(split, curve1, curve2, reachable1, reachable2, free_intervals1, free_intervals2);
            if (isLessThan) {
                ub = split;
            }
            else {
                lb = split;
            }
            if (Config::verbosity > 2) py::print("CFD: narrowed distance to to [", lb, ", ", ub, "]");
        }
    }
    
    const auto end = std::clock();
    result.value = ub;
    result.time_searches = (end - start) / CLOCKS_PER_SEC;
    result.number_searches = number_searches;
    return result;
}

bool _less_than_or_equal(const distance_t distance, Curve const& curve1, Curve const& curve2, 
        std::vector<Parameters> &reachable1, std::vector<Parameters> &reachable2,
        std::vector<Intervals> &free_intervals1, std::vector<Intervals> &free_intervals2) {
    
    if (Config::verbosity > 2) py::print("CFD: constructing FSD");
    const distance_t dist_sqr = distance * distance;
    const auto infty = std::numeric_limits<parameter_t>::infinity();
    const curve_size_t n1 = curve1.complexity();
    const curve_size_t n2 = curve2.complexity();

    if (Config::verbosity > 2) py::print("CFD: resetting old FSD");
    
    #pragma omp parallel for collapse(2) if (n1 * n2 > 1000)
    for (curve_size_t i = 0; i < n1; ++i) {
        for (curve_size_t j = 0; j < n2; ++j) {
            if (i < n1 - 1) reachable1[i][j] = infty;
            if (j < n2 - 1) reachable2[i][j] = infty;
            free_intervals1[j][i].reset();
            free_intervals2[i][j].reset();
        }
    }
    
    if (Config::verbosity > 2) py::print("CFD: FSD borders");
    
    for (curve_size_t i = 0; i < n1 - 1; ++i) {
        reachable1[i][0] = 0;
        if (curve2[0].dist_sqr(curve1[i+1]) > dist_sqr) break;
    }
    
    for (curve_size_t j = 0; j < n2 - 1; ++j) {
        reachable2[0][j] = 0;
        if (curve1[0].dist_sqr(curve2[j+1]) > dist_sqr) break;
    }
    
    if (Config::verbosity > 2) py::print("CFD: computing free space");
    
    #pragma omp parallel for collapse(2) if (n1 * n2 > 1000)
    for (curve_size_t i = 0; i < n1; ++i) {
        for (curve_size_t j = 0; j < n2; ++j) {
            if ((i < n1 - 1) and (j > 0)) {
                free_intervals1[j][i] = curve2[j].ball_intersection_interval(dist_sqr, curve1[i], curve1[i+1]);
            }
            if ((j < n2 - 1) and (i > 0)) {
                free_intervals2[i][j] = curve1[i].ball_intersection_interval(dist_sqr, curve2[j], curve2[j+1]);
            }
        }
    }
    
    if (Config::verbosity > 2) py::print("CFD: computing reachable space");
    
    for (curve_size_t i = 0; i < n1; ++i) {
        for (curve_size_t j = 0; j < n2; ++j) {
            if ((i < n1 - 1) and (j > 0)) {
                if (not free_intervals1[j][i].empty()) {
                    if (reachable2[i][j-1] != infty) {
                        reachable1[i][j] = free_intervals1[j][i].begin();
                    }
                    else if (reachable1[i][j-1] <= free_intervals1[j][i].end()) {
                        reachable1[i][j] = std::max(free_intervals1[j][i].begin(), reachable1[i][j-1]);
                    }
                }
            }
            if ((j < n2 - 1) and (i > 0)) {
                if (not free_intervals2[i][j].empty()) {
                    if (reachable1[i-1][j] != infty) {
                        reachable2[i][j] = free_intervals2[i][j].begin();
                    }
                    else if (reachable2[i-1][j] <= free_intervals2[i][j].end()) {
                        reachable2[i][j] = std::max(free_intervals2[i][j].begin(), reachable2[i-1][j]);
                    }
                }
            }
        }
    }
    return reachable1.back().back() < infty;
}

distance_t _greedy_upper_bound(const Curve &curve1, const Curve &curve2) {
    distance_t result = 0;
    
    const curve_size_t len1 = curve1.complexity(), len2 = curve2.complexity();
    curve_size_t i = 0, j = 0;
    
    while ((i < len1 - 1) and (j < len2 - 1)) {
        result = std::max(result, curve1[i].dist_sqr(curve2[j]));
        
        distance_t dist1 = curve1[i+1].dist_sqr(curve2[j]),
            dist2 = curve1[i].dist_sqr(curve2[j+1]),
            dist3 = curve1[i+1].dist_sqr(curve2[j+1]);
        
        if ((dist1 <= dist2) and (dist1 <= dist3)) ++i;
        else if ((dist2 <= dist1) and (dist2 <= dist3)) ++j;
        else {
            ++i;
            ++j;
        }
    }
    
    while (i < len1) result = std::max(result, curve1[i++].dist_sqr(curve2[j]));
    
    --i;
    
    while (j < len2) result = std::max(result, curve1[i].dist_sqr(curve2[j++]));
    
    return std::sqrt(result);
}

distance_t _projective_lower_bound(const Curve &curve1, const Curve &curve2) {
    std::vector<distance_t> distances1_sqr = std::vector<distance_t>(curve2.complexity() - 1), distances2_sqr = std::vector<distance_t>(curve1.complexity() + curve2.complexity() + 2);
    
    for (curve_size_t i = 0; i < curve1.complexity(); ++i) {
        #pragma omp parallel for
        for (curve_size_t j = 0; j < curve2.complexity() - 1; ++j) {
            if (curve2[j].dist_sqr(curve2[j+1]) > 0) {
                distances1_sqr[j] = curve1[i].line_segment_dist_sqr(curve2[j], curve2[j+1]);
            } else {
                distances1_sqr[j] = curve1[i].dist_sqr(curve2[j]);
            }
        }
        distances2_sqr[i] = *std::min_element(distances1_sqr.begin(), distances1_sqr.end());
    }
    
    distances1_sqr = std::vector<distance_t>(curve1.complexity() - 1);
    
    for (curve_size_t i = 0; i < curve2.complexity(); ++i) {
        #pragma omp parallel for
        for (curve_size_t j = 0; j < curve1.complexity() - 1; ++j) {
            if (curve1[j].dist_sqr(curve1[j+1]) > 0) {
                distances1_sqr[j] = curve2[i].line_segment_dist_sqr(curve1[j], curve1[j+1]);
            } else {
                distances1_sqr[j] = curve2[i].dist_sqr(curve1[j]);
            }
        }
        distances2_sqr[curve1.complexity() + i] = *std::min_element(distances1_sqr.begin(), distances1_sqr.end());
    }
    
    distances2_sqr[curve1.complexity() + curve2.complexity()] = curve1[0].dist_sqr(curve2[0]);
    distances2_sqr[curve1.complexity() + curve2.complexity() + 1] = curve1[curve1.complexity()-1].dist_sqr(curve2[curve2.complexity()-1]);
    return std::sqrt(*std::max_element(distances2_sqr.begin(), distances2_sqr.end()));
}

} // end namespace Continuous

namespace Discrete {
    
std::string Distance::repr() const {
    std::stringstream ss;
    ss << value;
    return ss.str();
}
    
Distance distance(const Curve &curve1, const Curve &curve2) {
    Distance result;
    const auto start = std::clock();
    
    std::vector<std::vector<distance_t>> a(curve1.complexity(), std::vector<distance_t>(curve2.complexity()));
    std::vector<std::vector<distance_t>> dists(curve1.complexity(), std::vector<distance_t>(curve2.complexity()));
    
    #pragma omp parallel for collapse(2)
    for (curve_size_t i = 0; i < curve1.complexity(); ++i) {
        for (curve_size_t j = 0; j < curve2.complexity(); ++j) {
            dists[i][j] = curve1[i].dist_sqr(curve2[j]);
        }
    }
    
    for (curve_size_t i = 0; i < curve1.complexity(); ++i) {
        for (curve_size_t j = 0; j < curve2.complexity(); ++j) {
            if (i == 0 and j == 0) a[i][j] = dists[i][j];
            else if (i == 0 and j > 0) a[i][j] = std::max(a[i][j-1], dists[i][j]);
            else if (i > 0 and j == 0) a[i][j] = std::max(a[i-1][j], dists[i][j]);
            else {
                a[i][j] = std::max(std::min(std::min(a[i-1][j], a[i-1][j-1]), a[i][j-1]), dists[i][j]);
            }
        }
    }
    
    const auto value = std::sqrt(a[curve1.complexity() - 1][curve2.complexity() - 1]);
    
    auto end = std::clock();
    
    result.time = (end - start) / CLOCKS_PER_SEC;
    result.value = value;
    return result;
    
}

} // end namespace Discrete

} // end namespace Frechet
