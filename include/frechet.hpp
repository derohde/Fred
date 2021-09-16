/*
Copyright 2020-2021 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include "config.hpp"
#include "types.hpp"
#include "point.hpp"
#include "interval.hpp"
#include "curve.hpp"

namespace Frechet {
namespace Continuous {
    
    extern distance_t error;
    
    struct Distance {
        distance_t value;
        double time_searches;
        double time_bounds;
        std::size_t number_searches;
        
        std::string repr() const;
    };
    
    Distance distance(const Curve&, const Curve&);
    
    Distance _distance(const Curve&, const Curve&, distance_t, distance_t);
            
    bool _less_than_or_equal(const distance_t, const Curve&, const Curve&, 
            std::vector<Parameters>&, std::vector<Parameters>&, 
            std::vector<Intervals>&, std::vector<Intervals>&);
            
    distance_t _greedy_upper_bound(const Curve&, const Curve&);
    distance_t _projective_lower_bound(const Curve&, const Curve&);
}
namespace Discrete {
    
    struct Distance {
        distance_t value;
        double time;
        
        std::string repr() const;
    };
    
    Distance distance(const Curve&, const Curve&);
    
    distance_t _dp(std::vector<Distances> &a, const curve_size_t i, const curve_size_t j, 
            const Curve &curve1, const Curve &curve2);
}
}
