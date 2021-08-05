 /*
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
 
#include "simplification.hpp"
 
Curve Simplification::approximate_weak_minimum_link_simplification(const Curve &pcurve, const distance_t epsilon) {
    Curve &curve = const_cast<Curve&>(pcurve);
    const curve_size_t complexity = curve.complexity();
    
    curve_size_t i = 0, j = 0;
    
    Curve simplification(curve.dimensions()), segment(2, curve.dimensions());
    simplification.push_back(curve.front());
    
    distance_t distance = 0;
    
    while (i < complexity - 1) {
        
        segment[0] = curve[i];
        j = 0;
                
        while (distance <= epsilon) {
            ++j;
            
            if (i + std::pow(2, j) >= complexity) break;
            
            curve.reset_subcurve();
            segment[1] = curve[i + std::pow(2, j)];
            curve.set_subcurve(i, i + std::pow(2, j));
            
            distance = Frechet::Continuous::distance(curve, segment).value;
        }
        
        curve_size_t low, mid, high;
        
        low = j == 1 ? 0 : std::pow(2, j - 1);
        high = std::min(static_cast<curve_size_t>(std::pow(2, j)), complexity - i - 1);
                
        while (low < high) {
            mid = std::ceil((low + high) / 2.);
                        
            curve.reset_subcurve();
            segment[1] = curve[i + mid];
            curve.set_subcurve(i, i + mid);
            
            distance = Frechet::Continuous::distance(curve, segment).value;
                                    
            if (distance <= epsilon) low = mid;
            else high = mid - 1;
        }
        
        i += low;
        curve.reset_subcurve();
        simplification.push_back(curve[i]);
    }
    return simplification;
}

Curve Simplification::approximate_weak_minimum_error_simplification(const Curve &curve, const curve_size_t ell) {
    Curve simplification(curve.dimensions()), segment(2, curve.dimensions());
    
    segment[0] = curve.front();
    segment[1] = curve.back();
    
    if (ell <= 2) return segment;
    
    distance_t min_distance = 0, max_distance = Frechet::Discrete::distance(curve, segment).value + 1, mid_distance;
    
    Curve new_simplification = Simplification::approximate_weak_minimum_link_simplification(curve, max_distance);

    while (new_simplification.complexity() > ell) {
        max_distance *= 2.;
        new_simplification = Simplification::approximate_weak_minimum_link_simplification(curve, max_distance);
    }
    
    while (max_distance - min_distance > Frechet::Continuous::epsilon) {
        mid_distance = (min_distance + max_distance) / 2.;
        
        new_simplification = Simplification::approximate_weak_minimum_link_simplification(curve, mid_distance);
        
        if (new_simplification.complexity() > ell) min_distance = mid_distance;
        else {
            simplification = new_simplification;
            max_distance = mid_distance;
        }
    }
    
    curve_size_t diff = ell - simplification.complexity();
    while (diff > 0) {
        simplification.push_back(simplification.back());
        --diff;
    }
    return simplification;
}
