/*
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#pragma once

#include "types.hpp"
#include "clustering.hpp"
#include "frechet.hpp"

namespace Coreset {
    
class Onemedian_Coreset {

    std::vector<curve_number_t> coreset;
    std::vector<distance_t> lambda;
    const distance_t Lambda = 76;
    distance_t cost;

public:
    Onemedian_Coreset() {}
    
    inline Onemedian_Coreset(const curve_size_t ell, const Curves &in, const distance_t epsilon, const double constant = 1) {
        compute(ell, in, epsilon, constant);
    }
    
    inline void compute(const curve_size_t ell, const Curves &in, const distance_t epsilon, const double eps, const bool round = true, const double constant = 1) {
        const auto n = in.size();
        const auto m = in.get_m();
        auto distances = Clustering::Distance_Matrix(in.size(), in.size());
        const auto c_approx = Clustering::arya(1, ell, in, distances, false);
        const auto center = c_approx.centers[0];
        cost = c_approx.value;
        if (cost == 0) {
            std::cerr << "WARNING: cost is zero, coreset construction not possible - check your input" << std::endl;
            return;
        }
        std::vector<double> probabilities(n);
        lambda = std::vector<distance_t>(n);
        
        for (curve_number_t i = 0; i < n; ++i) {
            lambda[i] = 52.0 / n + 24.0 / cost * Frechet::Continuous::distance(in[i], center).value;
            probabilities[i] = (lambda[i]) / Lambda;
        }
        
        auto prob_gen = Random::Custom_Probability_Generator<double>(probabilities);
        const std::size_t ssize = std::ceil(constant * 1/epsilon * 1/epsilon * std::log(m));
        const auto coreset_ind = prob_gen.get(ssize);
        for (curve_number_t i = 0; i < ssize; ++i) {
            coreset.push_back(coreset_ind[i]);
        }
    }
    
//     inline np::ndarray get_lambda() const {
//         np::dtype dt = np::dtype::get_builtin<distance_t>();
//         p::list l;
//         np::ndarray result = np::array(l, dt);
//         for (const auto &elem: lambda) {
//             l.append(elem);
//         }
//         result = np::array(l, dt);
//         return result;
//     }
    
    inline distance_t get_Lambda() const {
        return Lambda;
    }
    
//     inline np::ndarray get_curves() const {
//         np::dtype dt = np::dtype::get_builtin<curve_number_t>();
//         p::list l;
//         np::ndarray result = np::array(l, dt);
//         for (const auto &elem: coreset) {
//             l.append(elem);
//         }
//         result = np::array(l, dt);
//         return result;
//     }
    
    inline distance_t get_cost() const {
        return cost;
    }

};

};
