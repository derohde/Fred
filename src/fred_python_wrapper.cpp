/*
Copyright 2020-2021 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "config.hpp"
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

Curve minimum_error_simplification(const Curve &curve, const curve_size_t l) {
    Simplification::Subcurve_Shortcut_Graph graph(const_cast<Curve&>(curve));
    auto scurve = graph.minimum_error_simplification(l);
    scurve.set_name("Simplification of " + curve.get_name());
    return scurve;
}

Curve approximate_minimum_link_simplification(const Curve &curve, const distance_t epsilon) {
    auto scurve = Simplification::approximate_minimum_link_simplification(curve, epsilon);
    scurve.set_name("Simplification of " + curve.get_name());
    return scurve;
}

Curve approximate_minimum_error_simplification(const Curve &curve, const curve_size_t ell) {
    auto scurve = Simplification::approximate_minimum_error_simplification(curve, ell);
    scurve.set_name("Simplification of " + curve.get_name());
    return scurve;
}

PYBIND11_MODULE(backend, m) {
        
    py::class_<Config::Config>(m, "Config")
        .def(py::init<>())
        .def_property("continuous_frechet_error", [&](Config::Config&) { return fc::error; }, [&](Config::Config&, const distance_t error) { fc::error = error; })
        .def_property("verbosity", [&](Config::Config&) { return &Config::verbosity; }, [&](Config::Config&, const unsigned int verbosity) { Config::verbosity = verbosity; })
        .def_property("number_threads", [&](Config::Config&){ return &Config::number_threads; }, [&](Config::Config&, const int number_threads) {
            if (number_threads <= 0) {
                Config::number_threads = -1;
                Config::mp_dynamic = true;
            } else {
                Config::number_threads = number_threads;
                Config::mp_dynamic = false;
            }
#ifdef WITH_OMP
            omp_set_num_threads(Config::number_threads);
            omp_set_dynamic(Config::mp_dynamic);
#endif
        })
    ;
    
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
        .def("__add__", &Curves::operator+)
        .def("__len__", &Curves::number)
        .def("__str__", &Curves::str)
        .def("__iter__", [](Curves &v) { return py::make_iterator(v.begin(), v.end()); }, py::keep_alive<0, 1>())
        .def("__repr__", &Curves::repr)
        .def_property_readonly("values", &Curves::as_ndarray)
        .def_property_readonly("dimensions", &Curves::dimensions)
    ;
    
    py::class_<fc::Distance>(m, "Continuous_Frechet_Distance")
        .def_readwrite("time_searches", &fc::Distance::time_searches)
        .def_readwrite("time_bounds", &fc::Distance::time_bounds)
        .def_readwrite("number_searches", &fc::Distance::number_searches)
        .def_readwrite("value", &fc::Distance::value)
        .def("__repr__", &fc::Distance::repr)
    ;
    
    py::class_<fd::Distance>(m, "Discrete_Frechet_Distance")
        .def_readwrite("time", &fd::Distance::time)
        .def_readwrite("value", &fd::Distance::value)
        .def("__repr__", &fd::Distance::repr)
    ;
    
    py::class_<ddtw::Distance>(m, "Discrete_Dynamic_Time_Warping_Distance")
        .def_readwrite("time", &ddtw::Distance::time)
        .def_readwrite("value", &ddtw::Distance::value)
        .def("__repr__", &ddtw::Distance::repr)
    ;
    
    py::class_<Clustering::Clustering_Result>(m, "Clustering_Result")
        .def_readwrite("value", &Clustering::Clustering_Result::value)
        .def_readwrite("time", &Clustering::Clustering_Result::running_time)
        .def_readwrite("assignment", &Clustering::Clustering_Result::assignment)
        .def("__getitem__", &Clustering::Clustering_Result::get, py::return_value_policy::reference)
        .def("__setitem__", &Clustering::Clustering_Result::set)
        .def("__len__", &Clustering::Clustering_Result::size)
        .def("__iter__", [](Clustering::Clustering_Result &v) { return py::make_iterator(v.cbegin(), v.cend()); }, py::keep_alive<0, 1>())
        .def("compute_assignment", &Clustering::Clustering_Result::compute_assignment, py::arg("curves"), py::arg("consecutive_call") = false)
        .def("compute_center_enclosing_balls", &Clustering::Clustering_Result::compute_center_enclosing_balls)
    ;
    
    py::class_<Clustering::Cluster_Assignment>(m, "Cluster_Assignment")
        .def("__len__", &Clustering::Cluster_Assignment::size)
        .def("count", &Clustering::Cluster_Assignment::count)
        .def("get", &Clustering::Cluster_Assignment::get)
    ;
    
    py::class_<Coreset::Median_Coreset>(m, "Median_Coreset")
        .def(py::init<curve_number_t, curve_size_t, Curves&, parameter_t>())
        .def("cost", &Coreset::Median_Coreset::cost)
    ;
    
    m.def("continuous_frechet", &fc::distance);
    m.def("discrete_frechet", &fd::distance);
    m.def("discrete_dynamic_time_warping", &ddtw::distance);
    m.def("discrete_dynamic_time_warping_randomized", &ddtw::distance_randomized);
    
    m.def("minimum_error_simplification", &minimum_error_simplification);
    m.def("approximate_minimum_link_simplification", &approximate_minimum_link_simplification);
    m.def("approximate_minimum_error_simplification", &approximate_minimum_error_simplification);
    
    m.def("dimension_reduction", &JLTransform::transform_naive, py::arg("curves"), py::arg("epsilon") = 0.5, py::arg("empirical_constant") = true);

    m.def("discrete_klcenter", &Clustering::kl_center, py::arg("k") = 2, py::arg("l") = 2, py::arg("curves"), py::arg("local_search") = 0, py::arg("consecutive_call") = false, py::arg("random_start_center") = true, py::arg("fast_simplification") = false);
    m.def("discrete_klmedian", &Clustering::kl_median, py::arg("k") = 2, py::arg("l") = 2, py::arg("curves"), py::arg("consecutive_call") = false, py::arg("fast_simplification") = false);

}
