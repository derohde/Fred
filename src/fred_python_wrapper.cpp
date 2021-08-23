/*
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "curve.hpp"
#include "point.hpp"
#include "frechet.hpp"
#include "jl_transform.hpp"
#include "clustering.hpp"
#include "coreset.hpp"
#include "grid.hpp"
#include "simplification.hpp"
#include "dynamic_time_warping.hpp"

namespace py = pybind11;

namespace fc = Frechet::Continuous;
namespace fd = Frechet::Discrete;
namespace ddtw = Dynamic_Time_Warping::Discrete;

const distance_t default_epsilon = 0.001;

void set_frechet_epsilon(const double eps) {
    fc::epsilon = eps;
}

void set_frechet_rounding(const bool round) {
    fc::round = round;
}

distance_t get_frechet_epsilon() {
    return fc::epsilon;
}

bool get_frechet_rounding() {
    return fc::round;
}

Clustering::Clustering_Result klcenter(const curve_number_t num_centers, const curve_size_t ell, const Curves &in, Clustering::Distance_Matrix &distances, const Curves &center_domain = Curves(), const bool random_start_center = true) {
    auto result = Clustering::kl_center(num_centers, ell, in, distances, false, center_domain, random_start_center);
    return result;
}

Clustering::Clustering_Result klmedian(const curve_number_t num_centers, const curve_size_t ell, const Curves &in, Clustering::Distance_Matrix distances, const Curves &center_domain = Curves()) {

    auto result = Clustering::kl_median(num_centers, ell, in, distances, center_domain);
    
    return result;
}

// Clustering::Clustering_Result onemedian_sampling(const curve_size_t ell,  Curves &in, const double epsilon, const bool with_assignment = false, const Curves &center_domain = Curves()) {
//     
//     auto result = Clustering::one_median_sampling(ell, in, epsilon, with_assignment);
//     
//     return result;
// }
// 
// Clustering::Clustering_Result onemedian_exhaustive(const curve_size_t ell,  Curves &in, const bool with_assignment = false, const Curves &center_domain = Curves()) {
// 
//     auto result = Clustering::one_median_exhaustive(ell, in, with_assignment);
//     
//     return result;
// }
// 
// 
// Coreset::Onemedian_Coreset onemedian_coreset(const Curves &in, const curve_size_t ell, const double epsilon, const double constant = 1) {
//     return Coreset::Onemedian_Coreset(ell, in, epsilon, constant);
// }
// 

Curve weak_minimum_error_simplification(const Curve &curve, const curve_size_t l) {
    Simplification::Subcurve_Shortcut_Graph graph(const_cast<Curve&>(curve));
    auto scurve = graph.weak_minimum_error_simplification(l);
    scurve.set_name("Simplification of " + curve.get_name());
    return scurve;
}

Curve approximate_weak_minimum_link_simplification(const Curve &curve, const distance_t epsilon) {
    auto scurve = Simplification::approximate_weak_minimum_link_simplification(curve, epsilon);
    scurve.set_name("Simplification of " + curve.get_name());
    return scurve;
}

Curve approximate_weak_minimum_error_simplification(const Curve &curve, const curve_size_t ell) {
    auto scurve = Simplification::approximate_weak_minimum_error_simplification(curve, ell);
    scurve.set_name("Simplification of " + curve.get_name());
    return scurve;
}

void set_number_threads(std::uint64_t number) {
    omp_set_dynamic(0);
    omp_set_num_threads(number);
}

PYBIND11_MODULE(backend, m) {
    
    m.attr("default_epsilon_continuous_frechet") = default_epsilon;
    
    py::class_<Point>(m, "Point")
        .def(py::init<dimensions_t>())
        .def("__len__", &Point::dimensions)
        .def("__getitem__", &Point::get)
        .def("__setitem__", &Point::set)
        .def("__str__", &Point::str)
        .def("__iter__", [](Point &v) { return py::make_iterator(v.begin(), v.end()); }, py::keep_alive<0, 1>())
        .def("__repr__", &Point::repr)
        .def_property_readonly("values", &Point::as_ndarray)
    ;
    
    py::class_<Points>(m, "Points")
        .def(py::init<dimensions_t>())
        .def("__len__", &Points::number)
        .def("__getitem__", &Points::get, py::return_value_policy::reference)
        .def("__str__", &Points::str)
        .def("__iter__", [](Points &v) { return py::make_iterator(v.begin(), v.end()); }, py::keep_alive<0, 1>())
        .def("__repr__", &Points::repr)
        .def("add", &Points::add)
        .def_property_readonly("values", &Points::as_ndarray)
        .def_property_readonly("centroid", &Points::centroid)
    ;
    
    py::class_<Curve>(m, "Curve")
        .def(py::init<py::array_t<coordinate_t>>())
        .def(py::init<py::array_t<coordinate_t>, std::string>())
        .def_property_readonly("dimensions", &Curve::dimensions)
        .def_property_readonly("complexity", &Curve::complexity)
        .def_property("name", &Curve::get_name, &Curve::set_name)
        .def_property_readonly("values", &Curve::as_ndarray)
        .def_property_readonly("centroid", &Curve::centroid)
        .def("__getitem__", &Curve::get, py::return_value_policy::reference)
        .def("__len__", &Curve::complexity)
        .def("__str__", &Curve::str)
        .def("__iter__", [](Curve &v) { return py::make_iterator(v.begin(), v.end()); }, py::keep_alive<0, 1>())
        .def("__repr__", &Curve::repr)
    ;
        
    py::class_<Curves>(m, "Curves")
        .def(py::init<>())
        .def_property_readonly("m", &Curves::get_m)
        .def("add", &Curves::add)
        .def("simplify", &Curves::simplify)
        .def("__getitem__", &Curves::get, py::return_value_policy::reference)
        .def("__setitem__", &Curves::set)
        .def("__len__", &Curves::number)
        .def("__str__", &Curves::str)
        .def("__iter__", [](Curves &v) { return py::make_iterator(v.begin(), v.end()); }, py::keep_alive<0, 1>())
        .def("__repr__", &Curves::repr)
        .def_property_readonly("values", &Curves::as_ndarray)
    ;
    
    py::class_<fc::Distance>(m, "Continuous_Frechet_Distance")
        .def(py::init<>())
        .def_readwrite("time_searches", &fc::Distance::time_searches)
        .def_readwrite("time_bounds", &fc::Distance::time_bounds)
        .def_readwrite("number_searches", &fc::Distance::number_searches)
        .def_readwrite("value", &fc::Distance::value)
        .def("__repr__", &fc::Distance::repr)
    ;
    
    py::class_<fd::Distance>(m, "Discrete_Frechet_Distance")
        .def(py::init<>())
        .def_readwrite("time", &fd::Distance::time)
        .def_readwrite("value", &fd::Distance::value)
        .def("__repr__", &fd::Distance::repr)
    ;
    
    py::class_<ddtw::Distance>(m, "Discrete_Dynamic_Time_Warping_Distance")
        .def(py::init<>())
        .def_readwrite("time", &ddtw::Distance::time)
        .def_readwrite("value", &ddtw::Distance::value)
        .def("__repr__", &ddtw::Distance::repr)
    ;
    
    py::class_<Clustering::Distance_Matrix>(m, "Distance_Matrix")
        .def(py::init<>())
    ;
    
    py::class_<Clustering::Clustering_Result>(m, "Clustering_Result")
        .def(py::init<>())
        .def_readwrite("value", &Clustering::Clustering_Result::value)
        .def_readwrite("time", &Clustering::Clustering_Result::running_time)
        .def_readwrite("assignment", &Clustering::Clustering_Result::assignment)
        .def("__getitem__", &Clustering::Clustering_Result::get, py::return_value_policy::reference)
        .def("__len__", &Clustering::Clustering_Result::size)
        .def("__iter__", [](Clustering::Clustering_Result &v) { return py::make_iterator(v.cbegin(), v.cend()); }, py::keep_alive<0, 1>())
        .def("compute_assignment", &Clustering::Clustering_Result::compute_assignment)
    ;
    
    py::class_<Clustering::Cluster_Assignment>(m, "Cluster_Assignment")
        .def(py::init<>())
        .def("__len__", &Clustering::Cluster_Assignment::size)
        .def("count", &Clustering::Cluster_Assignment::count)
        .def("get", &Clustering::Cluster_Assignment::get)
    ;
    
    /*py::class_<Coreset::Onemedian_Coreset>("Onemedian_coreset")
        .def_property("lambd", &Coreset::Onemedian_Coreset::get_lambda)
        .def_property("Lambd", &Coreset::Onemedian_Coreset::get_Lambda)
        .def_property("cost", &Coreset::Onemedian_Coreset::get_cost)
        .def("curves", &Coreset::Onemedian_Coreset::get_curves)
    ;*/
    
    m.def("set_continuous_frechet_epsilon", &set_frechet_epsilon);
    m.def("set_continuous_frechet_rounding", &set_frechet_rounding);
    m.def("get_continuous_frechet_epsilon", &get_frechet_epsilon);
    m.def("get_continuous_frechet_rounding", &get_frechet_rounding);
    
    m.def("continuous_frechet", &fc::distance);
    m.def("discrete_frechet", &fd::distance);
    m.def("discrete_dynamic_time_warping", &ddtw::distance);
    
    m.def("weak_minimum_error_simplification", &weak_minimum_error_simplification);
    m.def("approximate_weak_minimum_link_simplification", &approximate_weak_minimum_link_simplification);
    m.def("approximate_weak_minimum_error_simplification", &approximate_weak_minimum_error_simplification);
    
    m.def("dimension_reduction", &JLTransform::transform_naive, py::arg("in") = Curves(), py::arg("epsilon") = 0.5, py::arg("empirical_constant") = true);

    m.def("discrete_klcenter", &klcenter, py::arg("num_centers") = 1, py::arg("ell") = 2, py::arg("in") = Curves(), py::arg("distances") = Clustering::Distance_Matrix(), py::arg("center_domain") = Curves(), py::arg("random_start_center") = true);
    m.def("discrete_klmedian", &klmedian, py::arg("num_centers") = 1, py::arg("ell") = 2, py::arg("in") = Curves(), py::arg("distances") = Clustering::Distance_Matrix(), py::arg("center_domain") = Curves());
    
    // these are experimental
    //m.def("two_two_dtw_one_two_median", &Clustering::two_two_dtw_one_two_median);
    //m.def("two_two_dtw_one_two_median_exact", &Clustering::two_two_dtw_one_two_median_exact);
    //def("discrete_onemedian_sampling", onemedian_sampling, onemedian_sampling_overloads());
    //def("discrete_onemedian_exhaustive", onemedian_exhaustive, onemedian_exhaustive_overloads());
    //def("onemedian_coreset", onemedian_coreset, onemedian_coreset_overloads());
    
    m.def("set_maximum_number_threads", set_number_threads);
}
