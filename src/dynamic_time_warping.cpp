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
    
Points vertices_matching_points(const Curve &curve1, const Curve &curve2, const Distance &dist) {
    if ((curve1.complexity() < 2) or (curve2.complexity() < 2)) {
        py::print("WARNING: curves must be of at least two points");
        Points result(curve1.dimensions());
        return result;
    }
    
    if (Config::verbosity > 1) py::print("DDTW: computing matching points from curve1 of complexity ", curve1.complexity(), " to curve2 of complexity ", curve2.complexity());
    if (Config::verbosity > 2) py::print("DDTW: distance between curve1 and curve2 is ", dist.value);
    
    std::vector<Points> matching_points(curve1.size(), Points(curve1.dimensions()));    
    
    for (curve_number_t i = 0; i < dist.matching.size(); ++i) {
        curve_number_t j = dist.matching[i].first, k = dist.matching[i].second;
        matching_points[j].push_back(curve2[k]);
    }
    
    Points result(curve1.size(), curve1.dimensions());
    
    for (curve_size_t i = 0; i < curve1.size(); ++i) {
        result[i] = matching_points[i].centroid();
    }
        
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
    std::vector<std::vector<distance_t>> a(curve1.complexity() + 1, std::vector<distance_t>(curve2.complexity() + 1, infty));
    std::vector<std::vector<distance_t>> dists(curve1.complexity(), std::vector<distance_t>(curve2.complexity()));
    
    #pragma omp parallel for collapse(2)
    for (curve_size_t i = 0; i < curve1.complexity(); ++i) {
        for (curve_size_t j = 0; j < curve2.complexity(); ++j) {
            dists[i][j] = curve1[i].dist(curve2[j]);
        }
    }
    
    unsigned int mwd = 0;
    
    curve_number_t ii = 1, jj = 1;
    distance_t min_ele;
    
    a[0][0] = 0;
    for (curve_size_t i = 1; i <= curve1.complexity(); ++i) {
        for (curve_size_t j = 1; j <= curve2.complexity(); ++j) {
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
    
    const auto end = std::clock();
    result.time = (end - start) / CLOCKS_PER_SEC;
    result.value = a[n1][n2];
        
    return result;
}

} // end namespace Discrete

} // end namespace Dynamic Time Warping
