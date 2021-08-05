/*
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "jl_transform.hpp"

namespace JLTransform {
    
Curves transform_naive(const Curves &in, const distance_t epsilon, const bool empirical_k = true) {
    
    if (in.empty()) return in;
    
    auto rg = Random::Gauss_Random_Generator<coordinate_t>(0, 1);
    
    curve_number_t number_points = 0;
    for (const Curve &elem : in) number_points += elem.complexity();
    
    const distance_t epsilonsq = epsilon * epsilon;       
    const distance_t epsiloncu = epsilonsq * epsilon;
    const dimensions_t new_number_dimensions = empirical_k ? std::ceil(2 * std::log(number_points) * 1/epsilonsq):
        std::ceil(4 * std::log(number_points) * 1 /((epsilonsq/2) - (epsiloncu/3)));
                                                
    std::vector<Coordinates> mat (new_number_dimensions);
    
    #if DEBUG
    std::cout << "populating " << new_number_dimensions << "x" << in[0].dimensions() << " matrix" << std::endl;
    #endif
    
    for (Coordinates &elem : mat) elem = rg.get(in[0].dimensions());
                
    Curves result(in.size(), in.get_m(), new_number_dimensions);
    
    coordinate_t sqrtk = std::sqrt(new_number_dimensions);
    
    #pragma omp parallel for
    for (curve_number_t l = 0; l < in.size(); ++l) {
        result[l] = Curve(in[l].complexity(), new_number_dimensions, in[l].get_name());
        
        for (curve_size_t i = 0; i < in[l].complexity(); ++i) {
            
            for (dimensions_t j = 0; j < new_number_dimensions; ++j) {
                result[l][i][j] = mat[j][0] * in[l][i][0];
                
                for (dimensions_t k = 1; k < in[l].dimensions(); ++k) {
                    result[l][i][j] += mat[j][k] * in[l][i][k];
                }
                
                result[l][i][j] /= sqrtk;
                
            }
            
        }
        
        #if DEBUG
        std::cout << "projected curve no. " << l << " from " << in.dimensions() << " to " << new_number_dimensions << " dimensions" << std::endl;
        #endif
        
    }
    
    return result;
}
    
}
