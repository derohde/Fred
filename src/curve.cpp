/*
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <typeinfo>

#include "curve.hpp"
#include "simplification.hpp"

Curve::Curve(const Points &points, const std::string &name) : Points(points), vstart{0}, vend{points.size() - 1}, name{name} {
    if (points.empty()) { 
        std::cerr << "warning: constructed empty curve" << std::endl;
        return; 
    }
    #if DEBUG
    std::cout << "constructed curve of complexity " << points.size() << std::endl;
    #endif
}

Curve::Curve(const py::array_t<coordinate_t> &in, const std::string &name) : Points(in.request().shape[0], in.request().ndim > 1 ? in.request().shape[1] : 1), name{name}, vstart{0}, vend{Points::size() - 1} {    
    const dimensions_t n_dimensions = in.ndim();
    auto shape = in.shape();
    const curve_size_t number_points = shape[0];
    
    if (n_dimensions > 2){
        std::cerr << "A Curve requires a 1- or 2-dimensional numpy array of type " << typeid(coordinate_t).name() << "." << std::endl;
        std::cerr << "Current dimensions: " << n_dimensions << std::endl;
        std::cerr << "WARNING: constructed empty curve" << std::endl;
        return;
    }

    if (n_dimensions == 2) {
        const dimensions_t point_size = shape[1];
        
        #if DEBUG
        std::cout << "constructing curve of size " << number_points << " and " << point_size << " dimensions" << std::endl;
        #endif
                
        #pragma omp parallel for simd
        for (curve_size_t i = 0; i < number_points; ++i) {
            for(dimensions_t j = 0; j < point_size; ++j){
              Points::operator[](i)[j] = *in.data(i, j);
            }
        }
    } else {
        #pragma omp parallel for simd
        for (curve_size_t i = 0; i < number_points; ++i) {
            Points::operator[](i)[0] = *in.data(i);
        }
    }
    
    if (empty()) { 
        std::cerr << "WARNING: constructed empty curve" << std::endl;
    return; 
    }
}

Curves Curves::simplify(const curve_size_t l) {
    const bool approx = true;
    Curves result(size(), l, Curves::dimensions());
    for (curve_number_t i = 0; i < size(); ++i) {
        if (approx) {
            Curve simplified_curve = Simplification::approximate_weak_minimum_error_simplification(std::vector<Curve>::operator[](i), l);
            simplified_curve.set_name("Simplification of " + std::vector<Curve>::operator[](i).get_name());
            result[i] = simplified_curve;
        } else {
            Simplification::Subcurve_Shortcut_Graph graph(std::vector<Curve>::operator[](i));
            Curve simplified_curve = graph.weak_minimum_error_simplification(l);
            simplified_curve.set_name("Simplification of " + std::vector<Curve>::operator[](i).get_name());
            result[i] = simplified_curve;
        }
        #if DEBUG
        std::cout << "Simplified curve " << i + 1 << "/" << size() << "." << std::endl;
        #endif
    }
    return result;
}

std::string Curve::repr() const {
    std::stringstream ss;
    ss << "fred.Curve '" << name << "' of complexity " << complexity() << " and " << dimensions() << " dimensions";
    return ss.str();
}

std::string Curves::repr() const {
    std::stringstream ss;
    ss << "fred.Curves collection with " << number() << " curves";
    return ss.str();
}

std::string Curve::str() const {
    std::stringstream ss;
    ss << name << std::endl;
    ss << *this;
    return ss.str();
}

std::string Curves::str() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

std::string Curve::get_name() const {
    return name;
}

void Curve::set_name(const std::string &name) {
    this->name = name;
}

std::ostream& operator<<(std::ostream &out, const Curve &curve) {
    if (curve.empty()) return out;
    out << "[";
    
    for (curve_size_t i = 0; i < curve.complexity() - 1; ++i) {
        out << curve[i] << ", ";
    }
    
    out << curve[curve.complexity() -1] << "]";

    return out;
}

std::ostream& operator<<(std::ostream &out, const Curves &curves) {
    if (curves.empty()) return out;
    out << "{";
    
    for (curve_number_t i = 0; i < curves.number() - 1; ++i) {
        out << curves[i] << ", ";
    }
    
    out << curves[curves.size() -1] << "}";

    return out;
}
