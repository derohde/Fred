/*
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <iostream> 

#include <boost/python/numpy.hpp>
#include <boost/python.hpp>

#include "point.hpp"
#include "interval.hpp"

namespace np = boost::python::numpy;
namespace p = boost::python;

class Curve {
public:
    typedef curve_size_t index_type;

    Curve(dimensions_t dimensions) : number_dimensions{dimensions} {}
    Curve(const Points& points, dimensions_t dimensions);
    Curve(const np::ndarray &in);

    inline std::size_t size() const { return points.size(); }
    inline bool empty() const { return points.empty(); }
    inline std::size_t dimensions() const { return number_dimensions; }
    inline Point const& operator[](const std::size_t i) const { return points[i]; }

    inline Point front() const { return points.front(); }
    inline Point back() const { return points.back(); }

    void push_back(const Point &point);

    inline Points::iterator begin() { return points.begin(); }
    inline Points::iterator end() { return points.end(); }
    inline Points::const_iterator begin() const { return points.cbegin(); }
    inline Points::const_iterator end() const { return points.cend(); }
    
private:
    const dimensions_t number_dimensions;
    Points points;
};

class Curves : public std::vector<Curve> {
public:
    inline Curve get(std::size_t i) const {
        return this->operator[](i);
    }
};

std::ostream& operator<<(std::ostream& out, const Curve& curve);
