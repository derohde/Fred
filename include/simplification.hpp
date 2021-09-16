/*
Copyright 2020-2021 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>

#include "config.hpp"
#include "types.hpp"
#include "curve.hpp"
#include "frechet.hpp"

namespace Simplification {

class Subcurve_Shortcut_Graph {

    Curve& curve;
    std::vector<Distances> edges;
    
public:
    
    Subcurve_Shortcut_Graph(Curve &curve) : curve{curve}, edges{std::vector<Distances>(curve.complexity(), Distances(curve.complexity(), std::numeric_limits<distance_t>::infinity()))} {
        if (Config::verbosity > 1) std::cout << "SIMPL: computing shortcut graph" << std::endl;
        const curve_size_t complexity = curve.complexity();
        Curve segment(2, curve.front().dimensions());
        auto distance = Frechet::Continuous::Distance();
        
        for (curve_size_t i = 0; i < complexity - 1; ++i) {
            for (curve_size_t j = i + 1; j < complexity; ++j) {
                
                if (Config::verbosity > 1) std::cout << "SIMPL: computing shortcut distance from vertex " << i << " to vertex " << j << std::endl;
                
                curve.set_subcurve(i, j);
                
                segment[0] = curve.front();
                segment[1] = curve.back();
                
                distance = Frechet::Continuous::distance(curve, segment);
                edges[i][j] = distance.value;
                
                curve.reset_subcurve();
            }
            
        }
    }
    
    Curve minimum_error_simplification(const curve_size_t ll) const {
        if (Config::verbosity > 1) std::cout << "SIMPL: computing exact minimum error simplification using shortcut graph" << std::endl;
        if (ll >= curve.complexity()) return curve;
        
        curve_size_t l = ll - 1;
        
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
                if (Config::verbosity > 1) std::cout << "SIMPL: initializing arrays" << std::endl;
                #pragma omp parallel for
                for (curve_size_t j = 1; j < curve.complexity(); ++j) {
                    distances[j][0] = edges[0][j];
                    predecessors[j][0] = 0;
                }
            } else {
                for (curve_size_t j = 1; j < curve.complexity(); ++j) {
                    if (Config::verbosity > 1) std::cout << "SIMPL: computing shortcut using " << i << " jumps" << std::endl;
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
        
        if (Config::verbosity > 1) std::cout << "SIMPL: backwards constructing simplification" << std::endl;
        
        curve_size_t ell = l;
        
        result.push_back(curve.back());
        curve_size_t predecessor = predecessors[curve.complexity() - 1][--ell];
        
        for (curve_size_t i = 0; i  < l; ++i) {
            result.push_back(curve[predecessor]);
            predecessor = predecessors[predecessor][--ell];
        }
                
        std::reverse(result.begin(), result.end());
        return result;
    }
    
};

Curve approximate_minimum_link_simplification(const Curve&, const distance_t);
Curve approximate_minimum_error_simplification(const Curve&, const curve_size_t);
 
};
