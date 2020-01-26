/*
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include "types.hpp"
#include "curve.hpp"
#include "random.hpp"

namespace JLTransform {
    
    Curves transform_naive(const Curves &in, const double epsilon, const bool empirical_k = true) {
        
        if (in.empty()) return in;
        
        auto rg = Gauss_Random_Generator<coordinate_t>(0, 1);
        
        auto number_points = 0;
        for (auto &elem: in) number_points += elem.size();
        
        const auto epsilonsq = epsilon * epsilon;       
        const auto epsiloncu = epsilonsq * epsilon;
        const auto new_number_dimensions = empirical_k ? std::ceil(2 * std::log(number_points) * 1/epsilonsq):
            std::ceil(4 * std::log(number_points) * 1 /((epsilonsq/2) - (epsiloncu/3)));
                                                    
        std::vector<std::vector<coordinate_t>> mat (new_number_dimensions);
        
        #if DEBUG
        std::cout << "populating " << new_number_dimensions << "x" << in[0].dimensions() << " matrix" << std::endl;
        #endif
        
        for (auto &elem: mat) elem = rg.get(in[0].dimensions());
                    
        Curves result;
        
        auto sqrtk = std::sqrt(new_number_dimensions);
        
        #pragma omp parallel for
        for (curve_size_t l = 0; l < in.size(); ++l) {
            result.push_back(Curve(new_number_dimensions));
            
            for (curve_size_t i = 0; i < in[l].size(); ++i) {
                std::vector<coordinate_t> new_coords(new_number_dimensions);
                
                for (curve_size_t j = 0; j < new_number_dimensions; ++j) {
                    new_coords[j] = mat[j][0] * in[l][i][0];
                    
                    for (curve_size_t k = 1; k < in[l].dimensions(); ++k) {
                        new_coords[j] += mat[j][k] * in[l][i][k];
                    }
                    
                    new_coords[j] /= sqrtk;
                    
                }
                
                Point new_point(new_coords);
                result[l].push_back(new_point);
            }
            
            #if DEBUG
            std::cout << "projected curve no. " << l << " from " << in[l].dimensions() << " to " << new_number_dimensions << " dimensions" << std::endl;
            #endif
            
        }
        
        return result;
    }
    
}
