/*
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <random>

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
    
    inline std::vector<T> get(const unsigned long n) {
        std::vector<T> results;
        for (auto i = n; i > 0; --i) {
            results.push_back(get());
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
    
    inline std::vector<T> get(const unsigned long n) {
        std::vector<T> results;
        for (auto i = n; i > 0; --i) {
            results.push_back(get());
        }
        return results;
    }
};

/*template <typename T = double>
class Rademacher_Random_Generator {
    std::mt19937_64 mersenne_twister;
    std::discrete_distribution<double> distribution;

public:
    Rademacher_Random_Generator() : mersenne_twister{std::random_device{}()}, distribution{-1/2, 1/2} {}
    
    inline double get() {
        return distribution(mersenne_twister);
    }
    
    inline std::vector<T> get(const unsigned long n) {
        std::vector<T> results;
        for (auto i = n; i > 0; --i) {
            results.push_back(get());
        }
        return results;
    }
};*/
