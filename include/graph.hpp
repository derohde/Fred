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
#include <map>
#include <algorithm>
#include <cmath>

#include "types.hpp"
#include "curve.hpp"
#include "frechet.hpp"

class Curve_Graph {

    Curve& curve;
    std::map<std::pair<curve_size_t, curve_size_t>, distance_t> edges;
    
public:
    
    Curve_Graph(Curve &curve) : curve{curve} {
        const auto complexity = curve.complexity();
        
        for (curve_size_t i = 0; i < complexity - 1; ++i) {
            for (curve_size_t j = i + 1; j < complexity; ++j) {
                curve.set_subcurve(i,j);
                Curve segment;
                segment.push_back(curve.front());
                segment.push_back(curve.back());
                auto distance = Frechet::Continuous::distance(curve, segment);
                edges.emplace(std::make_pair(i, j) , distance.value);
                curve.reset_subcurve();
            }
        }
    }
    
    distance_t cheapest_way_to(const curve_size_t l, const curve_size_t j) {
        if (j == 0 and l == 0) return 0;
        else if((j == 0 and l > 0) or l <= 0) return std::numeric_limits<distance_t>::infinity();
        else {
            std::vector<distance_t> andere(j);
            for (curve_size_t i = j-1; i > 0; --i) {
                auto search = edges.find(std::make_pair(i,j));
                andere[i] = search->second + cheapest_way_to(l - 1, i);
            }
            auto search = edges.find(std::make_pair(0,j));
            andere[0] = search->second;
            auto best = std::min_element(andere.begin(), andere.end());
            return andere[std::distance(andere.begin(), best)];
        }
    }
    
    Curve weak_minimum_error_simplification(const curve_size_t l) {
        Curve result;
        std::cout << cheapest_way_to(l-1, curve.complexity() - 1) << std::endl;
        return result;
    }
};
 
