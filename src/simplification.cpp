 /*
Copyright 2020-2023 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
 
#include "simplification.hpp"

namespace Frechet {

namespace Continuous {

namespace Simplification {
    
Subcurve_Shortcut_Graph::Subcurve_Shortcut_Graph(const Curve &pcurve) : curve{const_cast<Curve&>(pcurve)}, 
        edges{std::vector<std::vector<distance_t>>(curve.complexity(), std::vector<distance_t>(curve.complexity(), std::numeric_limits<distance_t>::infinity()))} {
            
    if (Config::verbosity > 1) py::print("SIMPL: computing shortcut graph");
    const curve_size_t complexity = curve.complexity();
    Curve segment(2, curve.front().dimensions());
    Distance dist;
    
    for (curve_size_t i = 0; i < complexity - 1; ++i) {
        for (curve_size_t j = i + 1; j < complexity; ++j) {
            
            if (Config::verbosity > 1) py::print("SIMPL: computing shortcut distance from vertex ", i, " to vertex ", j);
            
            curve.set_subcurve(i, j);
            
            segment[0] = curve.front();
            segment[1] = curve.back();
            
            dist = distance(curve, segment);
            edges[i][j] = dist.value;
            
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
    
    std::vector<std::vector<distance_t>> distances(curve.complexity(), std::vector<distance_t>(l, std::numeric_limits<distance_t>::infinity()));
    std::vector<std::vector<curve_size_t>> predecessors(curve.complexity(), std::vector<curve_size_t>(l));
    
    std::vector<distance_t> others;
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
    
    distance_t dist = 0;
    
    while (i < complexity - 1) {
        
        segment[0] = curve[i];
        j = 0;
        dist = 0;
        
        if (Config::verbosity > 1) py::print("ASIMPL: computing maximum length shortcut starting at ", i);
        
        if (Config::verbosity > 1) py::print("ASIMPL: exponential error search");
        
        while (dist <= epsilon) {
            ++j;
            
            if (i + std::pow(2, j) >= complexity) break;
            
            curve.reset_subcurve();
            segment[1] = curve[i + std::pow(2, j)];
            curve.set_subcurve(i, i + std::pow(2, j));
            
            dist = distance(curve, segment).value;
        }
        
        low = std::pow(2, j - 1);
        high = std::min(static_cast<curve_size_t>(std::pow(2, j)), complexity - i - 1);
        
        if (Config::verbosity > 1) py::print("ASIMPL: binary error search for low = ", low, " and high = ", high);
        
        while (low < high) {
            mid = std::ceil(low + (high - low) * .5);
                        
            curve.reset_subcurve();
            segment[1] = curve[i + mid];
            curve.set_subcurve(i, i + mid);
            
            dist = distance(curve, segment).value;
                                    
            if (dist > epsilon) high = mid - 1;
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
    if (ell >= curve.complexity()) return curve;
    Curve simplification(curve.dimensions()), segment(2, curve.dimensions());
    
    segment[0] = curve.front();
    segment[1] = curve.back();
    
    if (ell <= 2) return segment;
    
    distance_t min_distance = 0, max_distance = Frechet::Discrete::distance(curve, segment).value + 1, mid_distance;
    
    Curve new_simplification = approximate_minimum_link_simplification(curve, max_distance);

    if (Config::verbosity > 1) py::print("ASIMPL: computing upper bound for error by exponential search");
    while (new_simplification.complexity() > ell) {
        max_distance *= 2.;
        new_simplification = approximate_minimum_link_simplification(curve, max_distance);
    }
    
    if (Config::verbosity > 1) py::print("ASIMPL: binary search using upper bound");
    const distance_t epsilon = std::max(min_distance * Frechet::Continuous::error / 100, std::numeric_limits<distance_t>::epsilon());
    while (max_distance - min_distance > epsilon) {
        mid_distance = (min_distance + max_distance) / distance_t(2);
        if (mid_distance == max_distance or mid_distance == min_distance) break;
        new_simplification = approximate_minimum_link_simplification(curve, mid_distance);
        
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

}

}

}

namespace Dynamic_Time_Warping {

namespace Discrete {
    
namespace Simplification {
    
Curve approximate_minimum_error_simplification(const Curve &curve, const curve_size_t ell) {

    const auto dist_f = [&](const curve_size_t i, const curve_size_t j, const Point &x) {
        distance_t result = 0;
                        
        for (curve_size_t k = i; k <= j; ++k) {
            result += curve[k].dist(x);
        }
        
        return result;
    };
    
    const auto infty = std::numeric_limits<distance_t>::infinity();
    const curve_size_t n = curve.size(), m = n - 1;
    
    if (ell >= m) return curve;
    
    std::vector<std::vector<distance_t>> d(n, std::vector<distance_t>(ell, infty));
    std::vector<std::vector<Points>> c(n, std::vector<Points>(ell, Points(curve.dimensions())));
    
    distance_t tdist, mindist;
    curve_size_t minpoint, mip;
    
    // every prefix 0...i
    for (curve_size_t i = 0; i < curve.size(); ++i) {
        
        mindist = infty;
        
        // point that best summarizes prefix
        for (curve_size_t k = 0; k < curve.size(); ++k) {
            tdist = dist_f(0, i, curve[k]);
            
            if (tdist < mindist) {
                mindist = tdist;
                minpoint = k;
            }
        }
    
        c[i][0].push_back(curve[minpoint]);
        d[i][0] = mindist;
    }
        
    // for every infix
    for (curve_size_t i = 1; i < curve.size(); ++i) {
        
        // for every remaining simplification point
        for (curve_size_t j = 1; j < ell; ++j) {
            
        mindist = infty;
        
            // check for best prefix
            for (curve_size_t ip = 0; ip < i; ++ip) {
                
                // point that best summarizes prefix + infix
                for (curve_size_t k = 0; k < curve.size(); ++k) {
                    tdist = d[ip][j - 1] + dist_f(ip + 1, i, curve[k]);
                    
                    if (tdist < mindist) {
                        mindist = tdist;
                        minpoint = k;
                        mip = ip;
                    }
                }
                
            }
            
            c[i][j] = c[mip][j - 1];
            c[i][j].push_back(curve[minpoint]);
            d[i][j] = mindist;
                
        }
    }
    
    const curve_size_t min_ell = std::distance(d[m].begin(), std::min_element(d[m].begin(), d[m].end()));
    return Curve(c[m][min_ell]);
    
}

}
    
}
    
}
