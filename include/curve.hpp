/*
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <iostream> 
#include <string>
#include <sstream>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "types.hpp"
#include "point.hpp"
#include "interval.hpp"

namespace np = boost::python::numpy;
namespace p = boost::python;

class Curve : public Points {
public:    
    inline Curve() {}
    inline Curve(const curve_size_t m, const dimensions_t dimensions) : Points(m, Point(dimensions)) {}
    Curve(const Points &points);
    Curve(const np::ndarray &in);
    
    inline Point get(const curve_size_t i) const {
        return Points::operator[](i);
    }
    
    inline curve_size_t complexity() const { 
        return size(); 
    }
    
    inline dimensions_t dimensions() const { 
        return empty() ? 0 : Points::operator[](0).dimensions();
    }
    
    std::string str() const;
    
    std::string repr() const;
};

class Curves : public std::vector<Curve> {
    curve_size_t m;
    
public:
    Curves() {}
    Curves(const curve_number_t n, const curve_size_t m) : std::vector<Curve>(n), m{m} {}
    
    inline void add(Curve &curve) {
        push_back(curve);
        if (curve.complexity() > m) m = curve.complexity();
    }
    
    inline Curve get(const curve_number_t i) const {
        return std::vector<Curve>::operator[](i);
    }
    
    inline curve_size_t get_m() const {
        return m;
    }
    
    inline curve_number_t number() const {
        return size();
    }
    
    std::string str() const;
    
    std::string repr() const;
};

std::ostream& operator<<(std::ostream& out, const Curve&);
std::ostream& operator<<(std::ostream& out, const Curves&);
