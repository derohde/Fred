/*
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <random>
#include <algorithm>

namespace Random {

template <typename T = double>
class Uniform_Random_Generator {
    std::mt19937_64 mersenne_twister;
    std::uniform_real_distribution<T> distribution;
    T lbound, ubound;
    
public:
    Uniform_Random_Generator(const T lbound = 0, const T ubound = 1) : lbound{lbound}, ubound{ubound}, mersenne_twister{std::random_device{}()}, distribution{lbound, ubound} {}
    
    inline T get() {
        return distribution(mersenne_twister);
    }
    
    inline std::vector<T> get(const std::size_t n) {
        std::vector<T> results(n);
        for (std::size_t i = 0; i < n; ++i) {
            results[i] = get();
        }
        return results;
    }
};

template <typename T = double>
class Gauss_Random_Generator {
    std::mt19937_64 mersenne_twister;
    std::normal_distribution<T> distribution;
    
    T mean, stddev;
    
public:
    Gauss_Random_Generator(const T mean, const T stddev) : mean{mean}, stddev{stddev}, mersenne_twister{std::random_device{}()}, distribution{mean, stddev} {}
    
    inline T get() {
        return distribution(mersenne_twister);
    }
    
    inline std::vector<T> get(const std::size_t n) {
        std::vector<T> results(n);
        for (std::size_t i = 0; i < n; ++i) {
            results[i] = get();
        }
        return results;
    }
};

template <typename T = double>
class Custom_Probability_Generator {
    Uniform_Random_Generator<T> uform_gen;
    std::vector<T> cumulative_probabilities;
    
public:
    Custom_Probability_Generator(const std::vector<T> &probabilities) : uform_gen{0, 1} {
        if (probabilities.empty()) return;
        cumulative_probabilities = std::vector<T>(probabilities.size());
        cumulative_probabilities[0] = probabilities[0];
        for (std::size_t i = 1; i < probabilities.size(); ++i) 
            cumulative_probabilities[i] = probabilities[i] + cumulative_probabilities[i-1];
    }
    
    inline T get() {
        const std::size_t n = cumulative_probabilities.size();
        const T r = uform_gen.get();
        const auto upper = std::upper_bound(cumulative_probabilities.cbegin(), cumulative_probabilities.cend(), r);
        assert(upper != cumulative_probabilities.cend());
        const std::size_t result = std::distance(cumulative_probabilities.cbegin(), upper);
        return result;
    }
    
    inline std::vector<T> get(const std::size_t n) {
        std::vector<T> results(n);
        for (std::size_t i = 0; i < n; ++i) {
            results[i] = get();
        }
        return results;
    }
};

};
