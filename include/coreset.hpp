/*
Copyright 2020-2021 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#pragma once

#include "types.hpp"
#include "clustering.hpp"
#include "frechet.hpp"

namespace Coreset {
    
class Median_Coreset {
    const Curves &in;
    const curve_number_t k;
    const curve_size_t ell;
    const distance_t Lambda;
    parameter_t epsilon;
    distance_t constant;
    Clustering::Clustering_Result c_approx;
    Distances cluster_costs;
    Curve_Numbers cluster_sizes;
    Curve_Numbers coreset;
    Distances lambda;
    Parameters probabilities;

public:
    inline Median_Coreset(const curve_number_t k, curve_size_t ell, const Curves &in, const parameter_t epsilon, const distance_t constant = 1) : in{in}, k{k}, ell{ell}, epsilon{epsilon}, constant{constant}, cluster_costs(k, 0), cluster_sizes(k, 0), lambda(in.size()), Lambda{2*k + 12*std::sqrt(k) + 18}, probabilities(in.size()), c_approx{Clustering::kl_median(k, ell, in)} {
        c_approx.compute_assignment(in);
        for (curve_number_t i = 0; i < k; ++i) {
            for (curve_number_t j = 0; j < c_approx.assignment.count(i); ++j) {
                cluster_costs[i] += Frechet::Continuous::distance(in[c_approx.assignment.get(i, j)], c_approx.centers[i]).value;
                cluster_sizes[i]++;
            }
        }
        
        for (curve_number_t i = 0; i < k; ++i) {
            for (curve_number_t j = 0; j < c_approx.assignment.count(i); ++j) {
                lambda[c_approx.assignment.get(i, j)] = (1+ std::sqrt(2*k/18)) * (6 * Frechet::Continuous::distance(in[c_approx.assignment.get(i, j)], c_approx.centers[i]).value / c_approx.value + 6 * cluster_costs[i] / (c_approx.value * cluster_sizes[i])) + (1 + std::sqrt(18/(2*k))) * 2 / cluster_sizes[i];
                probabilities[c_approx.assignment.get(i, j)] = (lambda[c_approx.assignment.get(i, j)]) / Lambda;
            }
        }
        
        compute();
    }
    
    inline void compute() {
        const auto n = in.size();
        const auto m = in.get_m();
        const auto cost = c_approx.value;
        
        if (cost == 0) {
            std::cerr << "WARNING: cost is zero, coreset construction not possible - check your input" << std::endl;
            return;
        }
        
        auto prob_gen = Random::Custom_Probability_Generator<parameter_t>(probabilities);
        const std::size_t ssize = std::ceil(k * k * constant * 1/epsilon * 1/epsilon * std::log(m) * std::log(n));
        const auto coreset_ind = prob_gen.get(ssize);
        for (curve_number_t i = 0; i < ssize; ++i) {
            coreset.push_back(coreset_ind[i]);
        }
    }
    
    inline distance_t cost(const Curves &centers) const {
        distance_t result = 0;
        for (curve_size_t i = 0; i < coreset.size(); ++i) {
            distance_t min = std::numeric_limits<distance_t>::infinity();
            for (const auto &center : centers) {
                auto dist = Frechet::Continuous::distance(in[coreset[i]], center).value;
                if (dist < min) min = dist;
            }
            result += Lambda/(coreset.size() * lambda[coreset[i]]) * min;
        }
        return result;
    }

};

}
