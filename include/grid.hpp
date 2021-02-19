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
#include <cmath>

#include "types.hpp"
#include "point.hpp"

class Grid {
    
    Points P;
    
    Grid(Points P) : P{P} {};
    Grid(const Grid&);
    ~Grid() {};
    
public:
    
    static Grid build_cube_grid(const Point &p, const distance_t width, const distance_t edge_length) {
        const curve_size_t half_number_points_per_dimension = std::ceil(edge_length / width);
        const auto dimensions = p.dimensions();
        const curve_number_t number_points = std::pow(2 * half_number_points_per_dimension, dimensions);
        
        Coordinates coord1d(2 * half_number_points_per_dimension);
        for (curve_size_t i = 0; i < half_number_points_per_dimension; ++i) {
            coord1d[i] = -(half_number_points_per_dimension - i) * width;
            coord1d[half_number_points_per_dimension + i] = i * width; 
        }
        
        Points P(number_points, dimensions);
        auto counter = std::vector<curve_size_t>(dimensions, 0);
        
        for (curve_number_t i = 0; i < number_points; ++i) {
            P[i] = Point(dimensions);
            for (dimensions_t j = 0; j < dimensions; ++j) {
                
                P[i][j] = coord1d[counter[j]];
                
            }
            P[i] += p;
            for (dimensions_t j = 0; j < dimensions, i < P.size() - 1; ++j) {
                if (counter[j] < coord1d.size() - 1) {
                    ++counter[j];
                    break;
                } else {
                    counter[j] = 0;
                }
            }
        }
        
        return Grid{P};
    }
    
    const Points& get_points() const {
        return P;
    }
    
};
