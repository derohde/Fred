/*
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <boost/python.hpp>

#include "curve.hpp"
#include "point.hpp"
#include "frechet.hpp"
#include "jl_transform.hpp"
#include "clustering.hpp"
#include "coreset.hpp"
#include "grid.hpp"
#include "graph.hpp"

using namespace boost::python;
namespace np = boost::python::numpy;
namespace fc = Frechet::Continuous;
namespace fd = Frechet::Discrete;

const distance_t default_epsilon = 0.001;

fc::Result continuous_frechet(const Curve &curve1, const Curve &curve2) { 
    return fc::distance(curve1, curve2);
}

fd::Result discrete_frechet(const Curve &curve1, const Curve &curve2) {
    return fd::distance(curve1, curve2);
}

Curves jl_transform(const Curves &in, const double epsilon, const bool empirical_constant = true) {
    
    Curves curvesrp = JLTransform::transform_naive(in, epsilon, empirical_constant);
    
    return curvesrp;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(jl_transform_overloads, jl_transform, 2, 3);

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

Clustering::Clustering_Result klcenter(const curve_number_t num_centers, const curve_size_t ell, const Curves &in, const bool with_assignment = false) {
    auto result = Clustering::gonzalez(num_centers, ell, in, false, with_assignment);
    
    return result;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(klcenter_overloads, klcenter, 3, 4);

Clustering::Clustering_Result onemedian_sampling(const Curves &in, const double epsilon, const bool with_assignment = false) {
    
    auto result = Clustering::one_median_sampling(epsilon, in, with_assignment);
    
    return result;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(onemedian_sampling_overloads, onemedian_sampling, 2, 3);

Clustering::Clustering_Result onemedian_exhaustive(const Curves &in, const bool with_assignment = false) {

    auto result = Clustering::one_median_exhaustive(in, with_assignment);
    
    return result;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(onemedian_exhaustive_overloads, onemedian_exhaustive, 1, 2);

Clustering::Clustering_Result klmedian(const curve_number_t num_centers, const curve_size_t ell, const Curves &in, const bool with_assignment = false) {

    auto result = Clustering::arya(num_centers, ell, in, with_assignment);
    
    return result;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(klmedian_overloads, klmedian, 3, 4);

Coreset::Onemedian_Coreset onemedian_coreset(const Curves &in, const curve_size_t ell, const double epsilon, const double constant = 1) {
    return Coreset::Onemedian_Coreset(ell, in, epsilon, constant);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(onemedian_coreset_overloads, onemedian_coreset, 3, 4);

Curve weak_minimum_error_simplification(Curve &curve, curve_size_t l) {
    Subcurve_Graph graph(curve);
    return graph.weak_minimum_error_simplification(l);
}

BOOST_PYTHON_MODULE(backend)
{
    Py_Initialize();
    np::initialize();
    
    scope().attr("default_epsilon_continuous_frechet") = default_epsilon;
    
    class_<Point>("Point", init<>())
        .def("__len__", &Point::dimensions)
        .def("__getitem__", &Point::get)
        .def("__str__", &Point::str)
        .def("__iter__", range(&Point::cbegin, &Point::cend))
        .def("__repr__", &Point::repr)
        .add_property("values", &Point::as_ndarray)
    ;
    
    class_<Curve>("Curve", init<np::ndarray>())
        .def(init<np::ndarray, std::string>())
        .add_property("dimensions", &Curve::dimensions)
        .add_property("complexity", &Curve::complexity)
        .add_property("name", &Curve::get_name)
        .def("__getitem__", &Curve::get)
        .def("__len__", &Curve::complexity)
        .def("__str__", &Curve::str)
        .def("__iter__", range<return_value_policy<copy_const_reference>>(&Curve::cbegin, &Curve::cend))
        .def("__repr__", &Curve::repr)
    ;
    
    class_<Curves>("Curves", init<>())
        .add_property("m", &Curves::get_m)
        .def("add", &Curves::add)
        .def("__getitem__", &Curves::get)
        .def("__len__", &Curves::number)
        .def("__str__", &Curves::str)
        .def("__iter__", range<return_value_policy<copy_const_reference>>(&Curves::cbegin, &Curves::cend))
        .def("__repr__", &Curves::repr)
    ;
    
    class_<fc::Result>("Continuous_Frechet_Result", init<>())
        .add_property("time_searches", &fc::Result::time_searches)
        .add_property("time_bounds", &fc::Result::time_bounds)
        .add_property("number_searches", &fc::Result::number_searches)
        .add_property("value", &fc::Result::value)
        .def("__repr__", &fc::Result::repr)
    ;
    
    class_<fd::Result>("Discrete_Frechet_Result", init<>())
        .add_property("time", &fd::Result::time)
        .add_property("value", &fd::Result::value)
        .def("__repr__", &fd::Result::repr)
    ;
    
    class_<Clustering::Clustering_Result>("Clustering_Result", init<>())
        .add_property("value", &Clustering::Clustering_Result::value)
        .add_property("time", &Clustering::Clustering_Result::running_time)
        .add_property("assignment", &Clustering::Clustering_Result::assignment)
        .def("__getitem__", &Clustering::Clustering_Result::get)
        .def("__len__", &Clustering::Clustering_Result::size)
        .def("__iter__", range(&Clustering::Clustering_Result::cbegin, &Clustering::Clustering_Result::cend))
    ;
    
    class_<Clustering::Cluster_Assignment>("Clustering_Assignment", init<>())
        .def("__len__", &Clustering::Cluster_Assignment::size)
        .def("count", &Clustering::Cluster_Assignment::count)
        .def("get", &Clustering::Cluster_Assignment::get)
    ;
    
    class_<Coreset::Onemedian_Coreset>("Onemedian_coreset")
        .add_property("lambd", &Coreset::Onemedian_Coreset::get_lambda)
        .add_property("Lambd", &Coreset::Onemedian_Coreset::get_Lambda)
        .add_property("cost", &Coreset::Onemedian_Coreset::get_cost)
        .def("curves", &Coreset::Onemedian_Coreset::get_curves)
    ;
    
    def("set_continuous_frechet_epsilon", set_frechet_epsilon);
    def("set_continuous_frechet_rounding", set_frechet_rounding);
    def("get_continuous_frechet_epsilon", get_frechet_epsilon);
    def("get_continuous_frechet_rounding", get_frechet_rounding);
    def("dimension_reduction", jl_transform, jl_transform_overloads());
    def("continuous_frechet", continuous_frechet);
    def("discrete_frechet", discrete_frechet);
    def("discrete_klcenter", klcenter, klcenter_overloads());
    def("discrete_klmedian", klmedian, klmedian_overloads());
    def("discrete_onemedian_sampling", onemedian_sampling, onemedian_sampling_overloads());
    def("discrete_onemedian_exhaustive", onemedian_exhaustive, onemedian_exhaustive_overloads());
    def("onemedian_coreset", onemedian_coreset, onemedian_coreset_overloads());
    def("weak_minimum_error_simplification", weak_minimum_error_simplification);
}
