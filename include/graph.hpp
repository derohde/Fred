/*
Copyright 2020 Dennis Rohde

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

#include "types.hpp"
#include "curve.hpp"
#include "frechet.hpp"

class Subcurve_Graph {

    Curve& curve;
    std::vector<std::vector<distance_t>> edges;
    
public:
    
    Subcurve_Graph(Curve &curve) : curve{curve}, edges{std::vector<std::vector<distance_t>>(curve.complexity(), std::vector<distance_t>(curve.complexity(), std::numeric_limits<distance_t>::infinity()))} {
        const auto complexity = curve.complexity();
        
        for (curve_size_t i = 0; i < complexity - 1; ++i) {
            for (curve_size_t j = i + 1; j < complexity; ++j) {
                curve.set_subcurve(i, j);
                Curve segment;
                segment.push_back(curve.front());
                segment.push_back(curve.back());
                auto distance = Frechet::Continuous::distance(curve, segment);
                edges[i][j] = distance.value;
                curve.reset_subcurve();
            }
        }
    }
    
    Curve weak_minimum_error_simplification(const curve_size_t ll) {
        if (ll >= curve.complexity()) return curve;
        
        auto l = ll - 1;
        
        Curve result;
        
        if (ll <= 2) {
            result.push_back(curve.front());
            result.push_back(curve.back());
            return result;
        }
        
        std::vector<std::vector<distance_t>> distances(curve.complexity(), std::vector<distance_t>(l, std::numeric_limits<distance_t>::infinity()));
        std::vector<std::vector<curve_size_t>> predecessors(curve.complexity(), std::vector<curve_size_t>(l));
                
        for (curve_size_t i = 0; i < l; ++i) {
            
            if (i == 0) {
                
                for (curve_size_t j = 1; j < curve.complexity(); ++j) {
                    
                    distances[j][0] = edges[0][j];
                    predecessors[j][0] = 0;
                    
                }
                
            } else {
                
                for (curve_size_t j = 1; j < curve.complexity(); ++j) {
                    
                    std::vector<distance_t> others(j);
                    
                    for (curve_size_t k = 0; k < j; ++k) {
                        
                        others[k] = std::max(edges[k][j], distances[k][i - 1]);
                        
                    }
                    
                    auto best = std::distance(others.begin(), std::min_element(others.begin(), others.end()));
                    
                    distances[j][i] = others[best];
                    predecessors[j][i] = best;
                }
            }
        }
        
        auto ell = l;
        result.push_back(curve.back());
        auto predecessor = predecessors[curve.complexity() - 1][ell - 1];
        
        while(ell > 0) {
            result.push_back(curve[predecessor]);
            predecessor = predecessors[predecessor][--ell - 1];
        }
        
        std::reverse(result.begin(), result.end());
        return result;
    }
    
};
 
