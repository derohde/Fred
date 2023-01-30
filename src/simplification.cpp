 /*
Copyright 2020-2023 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
 
#include "simplification.hpp"

namespace Simplification {

Subcurve_Shortcut_Graph::Subcurve_Shortcut_Graph(const Curve &pcurve) : curve{const_cast<Curve&>(pcurve)}, 
        edges{std::vector<Distances>(curve.complexity(), Distances(curve.complexity(), std::numeric_limits<distance_t>::infinity()))} {
            
    if (Config::verbosity > 1) py::print("SIMPL: computing shortcut graph");
    const curve_size_t complexity = curve.complexity();
    Curve segment(2, curve.front().dimensions());
    auto distance = Frechet::Continuous::Distance();
    
    for (curve_size_t i = 0; i < complexity - 1; ++i) {
        for (curve_size_t j = i + 1; j < complexity; ++j) {
            
            if (Config::verbosity > 1) py::print("SIMPL: computing shortcut distance from vertex ", i, " to vertex ", j);
            
            curve.set_subcurve(i, j);
            
            segment[0] = curve.front();
            segment[1] = curve.back();
            
            distance = Frechet::Continuous::distance(curve, segment);
            edges[i][j] = distance.value;
            
            curve.reset_subcurve();
        }
        
    }
}

Curve Subcurve_Shortcut_Graph::minimum_error_simplification(const curve_size_t ll) const {
    if (Config::verbosity > 1) py::print("SIMPL: computing exact minimum error simplification using shortcut graph");
    if (ll >= curve.complexity()) return curve;
    
    const curve_size_t l = ll - 1;
    
    Curve result(curve.dimensions());
    
    if (ll <= 2) {
        result.push_back(curve.front());
        result.push_back(curve.back());
        return result;
    }
    
    std::vector<Distances> distances(curve.complexity(), Distances(l, std::numeric_limits<distance_t>::infinity()));
    std::vector<std::vector<curve_size_t>> predecessors(curve.complexity(), std::vector<curve_size_t>(l));
    
    Distances others;
    curve_size_t best = 0;
    for (curve_size_t i = 0; i < l; ++i) {
        
        if (i == 0) {
            if (Config::verbosity > 1) py::print("SIMPL: initializing arrays");
            #pragma omp parallel for
            for (curve_size_t j = 1; j < curve.complexity(); ++j) {
                distances[j][0] = edges[0][j];
                predecessors[j][0] = 0;
            }
        } else {
            for (curve_size_t j = 1; j < curve.complexity(); ++j) {
                if (Config::verbosity > 1) py::print("SIMPL: computing shortcut using ", i, " jumps");
                others.resize(j);
                #pragma omp parallel for
                for (curve_size_t k = 0; k < j; ++k) {
                    others[k] = std::max(distances[k][i - 1], edges[k][j]);
                }
                
                best = std::distance(others.begin(), std::min_element(others.begin(), others.end()));
                
                distances[j][i] = others[best];
                predecessors[j][i] = best;
            }
        }
    }
    
    if (Config::verbosity > 1) py::print("SIMPL: backwards constructing simplification");
    
    curve_size_t ell = l - 1;
    
    result.push_back(curve.back());
    curve_size_t predecessor = predecessors[curve.complexity() - 1][ell];
    
    for (curve_size_t i = 0; i < l - 1; ++i) {
        result.push_back(curve[predecessor]);
        predecessor = predecessors[predecessor][--ell];
    }
    
    result.push_back(curve.front());

    std::reverse(result.begin(), result.end());
    return result;
}
 
Curve approximate_minimum_link_simplification(const Curve &pcurve, const distance_t epsilon) {
    if (Config::verbosity > 1) py::print("ASIMPL: computing approximate minimum link simplification for curve of complexity ", pcurve.complexity());
    Curve &curve = const_cast<Curve&>(pcurve);
    const curve_size_t complexity = curve.complexity();
    
    curve_size_t i = 0, j = 0, low, mid, high;
    
    Curve simplification(curve.dimensions()), segment(2, curve.dimensions());
    simplification.push_back(curve.front());
    
    distance_t distance = 0;
    
    while (i < complexity - 1) {
        
        segment[0] = curve[i];
        j = 0;
        distance = 0;
        
        if (Config::verbosity > 1) py::print("ASIMPL: computing maximum length shortcut starting at ", i);
        
        if (Config::verbosity > 1) py::print("ASIMPL: exponential error search");
        
        while (distance <= epsilon) {
            ++j;
            
            if (i + std::pow(2, j) >= complexity) break;
            
            curve.reset_subcurve();
            segment[1] = curve[i + std::pow(2, j)];
            curve.set_subcurve(i, i + std::pow(2, j));
            
            distance = Frechet::Continuous::distance(curve, segment).value;
        }
        
        low = std::pow(2, j - 1);
        high = std::min(static_cast<curve_size_t>(std::pow(2, j)), complexity - i - 1);
        
        if (Config::verbosity > 1) py::print("ASIMPL: binary error search for low = ", low, " and high = ", high);
        
        while (low < high) {
            mid = std::ceil(low + (high - low) * .5);
                        
            curve.reset_subcurve();
            segment[1] = curve[i + mid];
            curve.set_subcurve(i, i + mid);
            
            distance = Frechet::Continuous::distance(curve, segment).value;
                                    
            if (distance > epsilon) high = mid - 1;
            else low = mid;
        }
        
        if (Config::verbosity > 1) py::print("ASIMPL: shortcutting from ", i, " to ", i+low);
        
        i += low;
        
        curve.reset_subcurve();
        
        simplification.push_back(curve[i]);
    }
    return simplification;
}

Curve approximate_minimum_error_simplification(const Curve &curve, const curve_size_t ell) {
    if (Config::verbosity > 1) py::print("ASIMPL: computing approximate minimum error simplification");
    Curve simplification(curve.dimensions()), segment(2, curve.dimensions());
    
    segment[0] = curve.front();
    segment[1] = curve.back();
    
    if (ell <= 2) return segment;
    
    distance_t min_distance = 0, max_distance = Frechet::Discrete::distance(curve, segment).value + 1, mid_distance;
    
    Curve new_simplification = Simplification::approximate_minimum_link_simplification(curve, max_distance);

    if (Config::verbosity > 1) py::print("ASIMPL: computing upper bound for error by exponential search");
    while (new_simplification.complexity() > ell) {
        max_distance *= 2.;
        new_simplification = Simplification::approximate_minimum_link_simplification(curve, max_distance);
    }
    
    if (Config::verbosity > 1) py::print("ASIMPL: binary search using upper bound");
    const distance_t epsilon = std::max(min_distance * Frechet::Continuous::error / 100, std::numeric_limits<distance_t>::epsilon());
    while (max_distance - min_distance > epsilon) {
        mid_distance = (min_distance + max_distance) / distance_t(2);
        if (mid_distance == max_distance or mid_distance == min_distance) break;
        new_simplification = Simplification::approximate_minimum_link_simplification(curve, mid_distance);
        
        if (new_simplification.complexity() > ell) min_distance = mid_distance;
        else {
            simplification = new_simplification;
            max_distance = mid_distance;
        }
    }
    if (Config::verbosity > 1) py::print("ASIMPL: backwards construction of simplification");
    curve_size_t diff = ell - simplification.complexity();
    while (diff > 0) {
        simplification.push_back(simplification.back());
        --diff;
    }
    return simplification;
}

};
