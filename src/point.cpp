/*
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "point.hpp"

std::string Point::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

std::string Point::repr() const {
    std::stringstream ss;
    ss << "fred.Point of " << dimensions() << " dimensions";
    return ss.str();
}

std::string Points::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

std::string Points::repr() const {
    std::stringstream ss;
    ss << size() << " fred.Points of " << dimensions() << " dimensions";
    return ss.str();
}

std::ostream& operator<<(std::ostream &out, const Point &p) {
    if (p.empty()) return out;
    out << "(";
    
    for (dimensions_t i = 0; i < p.dimensions() - 1; ++i){
        out << p[i] << ",";
    }
    
    out << p[p.dimensions() - 1] << ")";

    return out;
}

std::ostream& operator<<(std::ostream &out, const Points &p) {
    if (p.empty()) return out;
    out << "{";
    
    for (curve_size_t i = 0; i < p.size() - 1; ++i){
        out << p[i] << ",";
    }
    
    out << p[p.size() - 1] << "}";

    return out;
}
