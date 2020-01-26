/*
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <boost/chrono/include.hpp>

#include "types.hpp"
#include "point.hpp"
#include "interval.hpp"
#include "curve.hpp"

namespace Frechet {
namespace Continuous {
    
    struct Result {
        distance_t value;
        double time_searches;
        double time_bounds;
        std::size_t number_searches;
    };
    
    auto distance(const Curve&, const Curve&, const distance_t) -> Result;
    
    auto _greedyUpperBound(const Curve&, const Curve&) -> distance_t;
    
    auto _distance(const Curve&, const Curve&, const distance_t, const distance_t, 
            const distance_t = std::numeric_limits<distance_t>::epsilon()) -> Result;
            
    bool _lessThan(const distance_t, const Curve&, const Curve&, 
            std::vector<std::vector<distance_t>>&, std::vector<std::vector<distance_t>>&, 
            std::vector<std::vector<Interval>>&, std::vector<std::vector<Interval>>&);
}
namespace Discrete {
    
    struct Result {
        distance_t value;
        double time;
    };
    
    auto distance(const Curve&, const Curve&) -> Result;
    auto _dp(std::vector<std::vector<distance_t>> &a, const curve_size_t i, const curve_size_t j, 
            const Curve &curve1, const Curve &curve2) -> distance_t;
}
}
