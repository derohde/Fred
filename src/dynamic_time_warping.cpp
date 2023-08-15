/*
Copyright 2021 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <vector>
#include <limits>
#include <chrono>

#include "dynamic_time_warping.hpp"

namespace Dynamic_Time_Warping {

namespace Discrete {
        
std::string Distance::repr() const {
    std::stringstream ss;
    ss << value;
    return ss.str();
}
    
Points vertices_matching_points(const Curve &input_curve, const Curve &center_curve, const Distance &dist) {
    if ((input_curve.complexity() < 2) or (center_curve.complexity() < 2)) {
        py::print("WARNING: curves must be of at least two points");
        Points result(center_curve.dimensions());
        return result;
    }
        
    if (Config::verbosity > 1) py::print("DDTW: computing matching points from center_curve of complexity ", center_curve.complexity(), " to input_curve of complexity ", input_curve.complexity());
    if (Config::verbosity > 2) py::print("DDTW: distance between input_curve and center_curve is ", dist.value);
    
    std::vector<Points> matching_points(center_curve.size(), Points(center_curve.dimensions()));    
    
    curve_number_t j, k;
    
    for (curve_number_t i = 0; i < dist.matching.size(); ++i) {
        j = dist.matching[i].first;
        k = dist.matching[i].second;
        
        if (Config::verbosity > 2) py::print("DDTW: matching point ", j, " on input_curve to point ", k, " on curve 2");
        matching_points[k].push_back(input_curve[j]);
    }
    
    Points result(center_curve.size(), center_curve.dimensions());
    
    if (Config::verbosity > 2) py::print("DDTW: computing centroids to aggregate multi-matching points");
    
    for (curve_size_t i = 0; i < center_curve.size(); ++i) {
        if (Config::verbosity > 2) py::print("DDTW: computing centroid ", i);
        result[i] = matching_points[i].centroid();
    }

    if (Config::verbosity > 2) py::print("DDTW: matching points computed");
        
    return result;
}
    
Distance distance(const Curve &curve1, const Curve &curve2) {
    Distance result;
    
    if ((curve1.complexity() < 2) or (curve2.complexity() < 2)) {
        py::print("WARNING: curves must be of at least two points");
        return result;
    }
    
    const auto start = std::clock();
    
    const auto infty = std::numeric_limits<distance_t>::infinity();
        
    const curve_size_t n1 = curve1.complexity(), n2 = curve2.complexity();
    curve_size_t contingency1 = std::ceil(std::sqrt(n1)), contingency2 = std::ceil(std::sqrt(n2));
        
    if (n1 < n2) contingency1 += n2 - n1 + 1;
    if (n2 < n1) contingency2 += n1 - n2 + 1;

    std::vector<std::vector<std::pair<curve_number_t, curve_number_t>>> multi_warp_counter(n1 + 1, std::vector<std::pair<curve_number_t, curve_number_t>>(n2 + 1, std::make_pair(0, 0)));
    std::vector<std::vector<std::pair<curve_number_t, curve_number_t>>> b(n1 + 1, std::vector<std::pair<curve_number_t, curve_number_t>>(n2 + 1, std::make_pair(0, 0)));
    std::vector<std::vector<distance_t>> a(n1 + 1, std::vector<distance_t>(n2 + 1, infty));
    std::vector<std::vector<distance_t>> dists(n1, std::vector<distance_t>(n2));
    
    #pragma omp parallel for collapse(2)
    for (curve_size_t i = 0; i < n1; ++i) {
        for (curve_size_t j = 0; j < n2; ++j) {
            dists[i][j] = curve1[i].dist(curve2[j]);
        }
    }
    
    unsigned int mwd = 0;
    
    curve_number_t ii = 1, jj = 1;
    distance_t min_ele;
    
    a[0][0] = 0;
    for (curve_size_t i = 1; i <= n1; ++i) {
        for (curve_size_t j = 1; j <= n2; ++j) {
            a[i][j] = dists[i-1][j-1];
            mwd = 0;
            
            ii = i - 1, jj = j - 1;
            min_ele = a[i - 1][j - 1];
            
            if (a[i][j - 1] < min_ele) {
                if ((Config::dtw_contingency and multi_warp_counter[i][j - 1].first < contingency1) or not Config::dtw_contingency) {
                    min_ele = a[i][j - 1];
                    ii = i;
                    jj = j - 1;
                    mwd = 1;
                }
            }
            
            if (a[i - 1][j] < min_ele) {
                if ((Config::dtw_contingency and multi_warp_counter[i - 1][j].second < contingency2) or not Config::dtw_contingency) {
                    min_ele = a[i - 1][j];
                    ii = i - 1;
                    jj = j;
                    mwd = 2;
                }
            }
            
            multi_warp_counter[i][j] = multi_warp_counter[ii][jj];
            if (mwd == 1) ++multi_warp_counter[i][j].first;
            else if (mwd == 2) ++multi_warp_counter[i][j].second;
            
            a[i][j] += min_ele;
            b[i][j] = std::make_pair(ii, jj);
        }
    }
            
    for (curve_number_t i = n1, j = n2; ; i = b[i][j].first, j = b[i][j].second) {
        result.matching.push_back(std::make_pair(i - 1, j - 1));
        if (b[i][j].first == 1 and b[i][j].second == 1) break;
    }
    
    result.matching.push_back(std::make_pair(0, 0));
    result.n = n1;
    result.m = n2;
    
    const auto end = std::clock();
    result.time = (end - start) / CLOCKS_PER_SEC;
    result.value = a[n1][n2];
        
    return result;
}

} // end namespace Discrete

} // end namespace Dynamic Time Warping
