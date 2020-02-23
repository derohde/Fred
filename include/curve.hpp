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

class Curve {
    Points points;
    
public:
    typedef curve_size_t index_t;
    
    Curve() {}
    Curve(const curve_size_t m, const dimensions_t d) : points(m, Point(d)) {}
    Curve(const Points &points, const dimensions_t d);
    Curve(const np::ndarray &in);
    
    inline Point operator[](const index_t i) const {
        return points[i];
    }
    
    inline Point& operator[](const index_t i) {
        return points[i];
    }
    
    inline const Point& front() const { 
        return points.front();
    }
    
    inline const Point& back() const {
        return points.back();
    }
    
    inline Points::iterator begin() { 
        return points.begin(); 
    }
    
    inline Points::iterator end() { 
        return points.end(); 
    }
    
    inline Points::const_iterator cbegin() const {
        return points.cbegin();
    }
    
    inline Points::const_iterator cend() const {
        return points.cend();
    }
    
    inline curve_size_t size() const { 
        return points.size(); 
    }
    
    inline bool empty() const { 
        return points.empty(); 
    }
    
    inline dimensions_t dimensions() const { 
        if (points.empty()) return 0;
        else return points[0].size();
    }
    
    std::string str() const;
};

class Curves : public std::vector<Curve> {
    curve_size_t m;
    
public:
    Curves() {}
    Curves(const curve_size_t n, const curve_size_t m) : std::vector<Curve>{n}, m{m} {}
    
    inline void add(Curve &curve) {
        push_back(curve);
        if (curve.size() > m) m = curve.size();
    }
    
    inline Curve operator[](const curve_size_t i) const {
        return std::vector<Curve>::operator[](i);
    }
    
    inline curve_size_t get_m() const {
        return m;
    }
    
    inline curve_size_t number() const {
        return size();
    }
    
    virtual std::vector<Curve>::const_iterator cbegin() const {
        return std::vector<Curve>::cbegin();
    }
    
    virtual std::vector<Curve>::const_iterator cend() const {
        return std::vector<Curve>::cend();
    }
    
    std::string str() const;
};

std::ostream& operator<<(std::ostream& out, const Curve&);
std::ostream& operator<<(std::ostream& out, const Curves&);
