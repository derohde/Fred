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

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "types.hpp"
#include "point.hpp"
#include "interval.hpp"

namespace py = pybind11;

class Curve : private Points {
    
    curve_size_t vstart = 0, vend = 0;
    std::string name;
    
public:
    typedef typename Points::iterator iterator;
    
    inline Curve(const dimensions_t dim, const std::string &name = "unnamed curve") : Points(dim), vstart{0}, vend{0}, name{name} {}
    inline Curve(const curve_size_t m, const dimensions_t dimensions, const std::string &name = "unnamed curve") : Points(m, Point(dimensions)), vstart{0}, vend{m-1}, name{name} {}
    Curve(const Points &points, const std::string &name = "unnamed curve");
    Curve(const py::array_t<coordinate_t> &in, const std::string &name = "unnamed curve");
    
    inline Point& get(const curve_size_t i) {
        return Points::operator[](vstart + i);
    }
    
    inline const Point& operator[](const curve_size_t i) const {
        return Points::operator[](vstart + i);
    }
    
    inline Point& operator[](const curve_size_t i) {
        return Points::operator[](vstart + i);
    }
    
    inline const Point& front() const {
        return Points::operator[](vstart);
    }
    
    inline Point& front() {
        return Points::operator[](vstart);
    }
    
    inline const Point& back() const {
        return Points::operator[](vend);
    }
    
    inline Point& back() {
        return Points::operator[](vend);
    }
    
    inline const auto begin() const {
        return Points::begin();
    }
    
    inline const auto end() const {
        return Points::end();
    }
    
    inline auto begin() {
        return Points::begin();
    }
    
    inline auto end() {
        return Points::end();
    }
    
    inline bool empty() const {
        return Points::empty();
    }
    
    inline curve_size_t complexity() const {
        return empty() ? 0 : vend - vstart + 1; 
    }
    
    inline curve_size_t size() const {
        return empty() ? 0 : vend - vstart + 1;
    }
    
    inline dimensions_t dimensions() const { 
        return empty() ? 0 : Points::dimensions();
    }
    
    inline void set_subcurve(const curve_size_t start, const curve_size_t end) {
        vstart = start;
        vend = end;
    }
    
    inline void reset_subcurve() {
        vstart = 0;
        vend = Points::size() - 1;
    }
    
    inline void push_back(const Point &point) {
        Points::push_back(point);
        vend = Points::size() - 1;
    }
    
    inline void push_back(Point &&point) {
        Points::push_back(point);
        vend = Points::size() - 1;
    }
    
    inline Point centroid() const {
        return Points::centroid();
    }
    
    inline auto as_ndarray() const {
        py::list l;
        for (const Point &elem : *this) {
            l.append(elem.as_ndarray());
        }
        return py::array_t<coordinate_t>(l);
    }
    
    void set_name(const std::string&);
    
    std::string get_name() const;
    
    std::string str() const;
    
    std::string repr() const;
};

class Curves : public std::vector<Curve> {
    curve_size_t m;
    dimensions_t dim;
    
public:
    Curves(const dimensions_t dim = 0) : dim{dim} {}
    Curves(const curve_number_t n, const curve_size_t m, const dimensions_t dim) : std::vector<Curve>(n, Curve(dim)), m{m}, dim{dim} {}
    
    inline void add(Curve &curve) {
        if (curve.dimensions() != dim) {
            if (dim == 0) {
                dim = curve.dimensions();
            } else {
                std::cerr << "Wrong number of dimensions; expected " << dim << " dimensions and got " << curve.dimensions() << " dimensions." << std::endl;
                return;
            }
        }
        push_back(curve);
        if (curve.complexity() > m) m = curve.complexity();
    }
    
    inline Curve& get(const curve_number_t i) {
        return std::vector<Curve>::operator[](i);
    }
    
    inline void set(const curve_number_t i, const Curve &val) {
        std::vector<Curve>::operator[](i) = val;
    }
    
    inline curve_size_t get_m() const {
        return m;
    }
    
    inline curve_number_t number() const {
        return size();
    }
    
    inline dimensions_t dimensions() const {
        return dim;
    }
    
    inline auto as_ndarray() const {
        py::list l;
        for (const Curve &elem : *this) {
            l.append(elem.as_ndarray());
        }
        return py::array_t<coordinate_t>(l);
    }
    
    Curves simplify(const curve_size_t);
    
    std::string str() const;
    
    std::string repr() const;
};

std::ostream& operator<<(std::ostream& out, const Curve&);
std::ostream& operator<<(std::ostream& out, const Curves&);
