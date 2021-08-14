/*
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <vector>
#include <limits>
#include <chrono>

#include "frechet.hpp"

namespace Frechet {

namespace Continuous {
    
distance_t epsilon = 0.001;
bool round = true;
    
std::string Distance::repr() const {
    std::stringstream ss;
    ss << value;
    return ss.str();
}

Distance distance(const Curve &curve1, const Curve &curve2) {
    if ((curve1.complexity() < 2) or (curve2.complexity() < 2)) {
        std::cerr << "WARNING: comparison possible only for curves of at least two points" << std::endl;
        Distance result;
        result.value = std::numeric_limits<distance_t>::signaling_NaN();
        return result;
    }
    if (curve1.dimensions() != curve2.dimensions()) {
        std::cerr << "WARNING: comparison possible only for curves of equal number of dimensions" << std::endl;
        Distance result;
        result.value = std::numeric_limits<distance_t>::signaling_NaN();
        return result;
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    const distance_t lb = std::sqrt(std::max(curve1[0].dist_sqr(curve2[0]), curve1[curve1.complexity()-1].dist_sqr(curve2[curve2.complexity()-1])));
    const distance_t ub = _greedy_upper_bound(curve1, curve2);
    auto end = std::chrono::high_resolution_clock::now();

    #if DEBUG
    std::cout << "narrowed to [" << lb << ", " << ub.value << "]" << std::endl;
    #endif

    auto dist = _distance(curve1, curve2, ub, lb);
    dist.time_bounds = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    if (round) dist.value =  std::round(dist.value * 1e3) / 1e3;

    return dist;
}

Distance _distance(const Curve &curve1, const Curve &curve2, distance_t ub, distance_t lb) {
    Distance result;
    auto start = std::chrono::high_resolution_clock::now();
    
    distance_t split = (ub + lb)/2;
    std::size_t number_searches = 0;
    
    if (ub - lb > epsilon) {
        auto infty = std::numeric_limits<parameter_t>::infinity();
        std::vector<std::vector<parameter_t>> reachable1(curve1.complexity()-1, std::vector<parameter_t>(curve2.complexity(), infty));
        std::vector<std::vector<parameter_t>> reachable2(curve1.complexity(), std::vector<parameter_t>(curve2.complexity()-1, infty));
        
        std::vector<std::vector<Interval>> free_intervals1(curve2.complexity(), std::vector<Interval>(curve1.complexity(), Interval()));
        std::vector<std::vector<Interval>> free_intervals2(curve1.complexity(), std::vector<Interval>(curve2.complexity(), Interval()));

        if (std::isnan(lb) or std::isnan(ub)) {
            result.value = std::numeric_limits<distance_t>::signaling_NaN();
            return result;
        }

        //Binary search over the feasible distances
        while (ub - lb > epsilon) {
            ++number_searches;
            split = (ub + lb)/2;
            auto isLessThan = _less_than_or_equal(split, curve1, curve2, reachable1, reachable2, free_intervals1, free_intervals2);
            if (isLessThan) {
                ub = split;
            }
            else {
                lb = split;
            }
            #if DEBUG
            std::cout << "narrowed to [" << lb << ", " << ub << "]" << std::endl;
            #endif
        }
    }
    
    distance_t value = (ub + lb)/2.;
    auto end = std::chrono::high_resolution_clock::now();
    result.value = value;
    result.time_searches = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    result.number_searches = number_searches;
    return result;
}

bool _less_than_or_equal(const distance_t distance, Curve const& curve1, Curve const& curve2, 
        std::vector<std::vector<parameter_t>> &reachable1, std::vector<std::vector<parameter_t>> &reachable2,
        std::vector<std::vector<Interval>> &free_intervals1, std::vector<std::vector<Interval>> &free_intervals2) {
    
    const distance_t dist_sqr = distance * distance;
    auto infty = std::numeric_limits<parameter_t>::infinity();

    for (auto &elem: reachable1) {
        #pragma omp parallel for
        for (curve_size_t i = 0; i < elem.size(); ++i) {
            elem[i] = infty;
        }
    }
    
    for (auto &elem: reachable2) {
        #pragma omp parallel for
        for (curve_size_t i = 0; i < elem.size(); ++i) {
            elem[i] = infty;
        }
    }
    
    for (auto &elem: free_intervals1) {
        #pragma omp parallel for
        for (curve_size_t i = 0; i < elem.size(); ++i) {
            elem[i].reset();
        }
    }
    
    for (auto &elem: free_intervals2) {
        #pragma omp parallel for
        for (curve_size_t i = 0; i < elem.size(); ++i) {
            elem[i].reset();
        }
    }
    
    for (curve_size_t i = 0; i < curve1.complexity() - 1; ++i) {
        reachable1[i][0] = 0;
        if (curve2[0].dist_sqr(curve1[i+1]) > dist_sqr) break;
    }
    
    for (curve_size_t j = 0; j < curve2.complexity() - 1; ++j) {
        reachable2[0][j] = 0;
        if (curve1[0].dist_sqr(curve2[j+1]) > dist_sqr) break;
    }
    
    #pragma omp parallel for collapse(2)
    for (curve_size_t i = 0; i < curve1.complexity(); ++i) {
        for (curve_size_t j = 0; j < curve2.complexity(); ++j) {
            if ((i < curve1.complexity() - 1) and (j > 0)) {
                free_intervals1[j][i] = curve2[j].intersection_interval(dist_sqr, curve1[i], curve1[i+1]);
            }
            if ((j < curve2.complexity() - 1) and (i > 0)) {
                free_intervals2[i][j] = curve1[i].intersection_interval(dist_sqr, curve2[j], curve2[j+1]);
            }
        }
    }
    
    for (curve_size_t i = 0; i < curve1.complexity(); ++i) {
        for (curve_size_t j = 0; j < curve2.complexity(); ++j) {
            if ((i < curve1.complexity() - 1) and (j > 0)) {
                if (not free_intervals1[j][i].is_empty()) {
                    if (reachable2[i][j-1] != infty) {
                        reachable1[i][j] = free_intervals1[j][i].begin();
                    }
                    else if (reachable1[i][j-1] <= free_intervals1[j][i].end()) {
                        reachable1[i][j] = std::max(free_intervals1[j][i].begin(), reachable1[i][j-1]);
                    }
                }
            }
            if ((j < curve2.complexity() - 1) and (i > 0)) {
                if (not free_intervals2[i][j].is_empty()) {
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

} // end namespace Continuous

namespace Discrete {
    
std::string Distance::repr() const {
    std::stringstream ss;
    ss << value;
    return ss.str();
}
    
Distance distance(const Curve &curve1, const Curve &curve2) {
    Distance result;
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<distance_t>> a(curve1.complexity(), std::vector<distance_t>(curve2.complexity(), -1));
    auto value = std::sqrt(_dp(a, curve1.complexity() - 1, curve2.complexity() - 1, curve1, curve2));
    auto end = std::chrono::high_resolution_clock::now();
    result.time = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    result.value = value;
    return result;
    
}

distance_t _dp(std::vector<std::vector<distance_t>> &a, const curve_size_t i, const curve_size_t j, const Curve &curve1, const Curve &curve2) {
    if (a[i][j] > -1) return a[i][j];
    else if (i == 0 and j == 0) return curve1[i].dist_sqr(curve2[j]);
    else if (i > 0 and j == 0) return std::max(_dp(a, i-1, 0, curve1, curve2), curve1[i].dist_sqr(curve2[j]));
    else if (i == 0 and j > 0) return std::max(_dp(a, 0, j-1, curve1, curve2), curve1[i].dist_sqr(curve2[j]));
    else {
        a[i][j] = std::max(
                    std::min(
                        std::min(_dp(a, i-1, j, curve1, curve2), 
                            _dp(a, i-1, j-1, curve1, curve2)), 
                        _dp(a, i, j-1, curve1, curve2)), 
                    curve1[i].dist_sqr(curve2[j]));
    }
    return a[i][j];
}

} // end namespace Discrete

} // end namespace Frechet
