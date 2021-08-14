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
    
Distance distance(const Curve &curve1, const Curve &curve2) {
    Distance result;
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<distance_t>> a(curve1.complexity() + 1, std::vector<distance_t>(curve2.complexity() + 1, std::numeric_limits<distance_t>::infinity()));
    a[0][0] = 0;
    for (curve_size_t i = 1; i <= curve1.complexity(); ++i) {
        for (curve_size_t j = 1; j <= curve2.complexity(); ++j) {
            a[i][j] = std::sqrt(curve1[i-1].dist_sqr(curve2[j-1]));
            a[i][j] += std::min(std::min(a[i-1][j], a[i][j-1]), a[i-1][j-1]);
        }
    }
    auto value = a[curve1.complexity()][curve2.complexity()];
    auto end = std::chrono::high_resolution_clock::now();
    result.time = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    result.value = value;
    return result;
    
}

} // end namespace Discrete

} // end namespace Dynamic Time Warping
