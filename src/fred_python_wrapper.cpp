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

using namespace boost::python;
namespace np = boost::python::numpy;
namespace fc = Frechet::Continuous;
namespace fd = Frechet::Discrete;

distance_t epss = 0.001;

fc::Result continuous_frechet(const Curve &curve1, const Curve &curve2, const distance_t eps = 0.001, const bool round = true) { 
    return fc::distance(curve1, curve2, eps, round);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(continuous_frechet_overloads, continuous_frechet, 2, 4);

fd::Result discrete_frechet(const Curve &curve1, const Curve &curve2) {
    return fd::distance(curve1, curve2);
}

Curves jl_transform(const Curves &in, const double epsilon, const bool empirical_constant = true) {
    
    Curves curvesrp = JLTransform::transform_naive(in, epsilon, empirical_constant);
    
    return curvesrp;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(jl_transform_overloads, jl_transform, 2, 3);

Clustering::Clustering_Result kcenter(const curve_size_t num_centers, const Curves &in, 
                                        const distance_t eps = 0.001, const bool round = true, const bool with_assignment = false) {
    auto result = Clustering::gonzalez(num_centers, in, eps, round, false, with_assignment);
    
    return result;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(kcenter_overloads, kcenter, 2, 5);

Clustering::Clustering_Result onemedian_sampling(const Curves &in, const double epsilon, 
                                                    const distance_t eps = 0.001, const bool round = true, const bool with_assignment = false) {
    
    auto result = Clustering::one_median_sampling(epsilon, in, eps, round, with_assignment);
    
    return result;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(onemedian_sampling_overloads, onemedian_sampling, 2, 5);

Clustering::Clustering_Result onemedian_exhaustive(const Curves &in, 
                                                    const distance_t eps = 0.001, const bool round = true, const bool with_assignment = false) {

    auto result = Clustering::one_median_exhaustive(in, eps, round, with_assignment);
    
    return result;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(onemedian_exhaustive_overloads, onemedian_exhaustive, 1, 4);

Clustering::Clustering_Result kmedian(const curve_size_t num_centers, const Curves &in, 
                                        const distance_t eps = 0.001, const bool round = true, const bool with_assignment = false) {

    auto result = Clustering::arya(num_centers, in, eps, round, with_assignment);
    
    return result;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(kmedian_overloads, kmedian, 2, 5);

Coreset::Onemedian_Coreset onemedian_coreset(const Curves &in, const double epsilon, 
                                                const double eps = 0.001, const bool round = true, const double constant = 1) {
    return Coreset::Onemedian_Coreset(in, epsilon, eps, round, constant);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(onemedian_coreset_overloads, onemedian_coreset, 2, 5);

BOOST_PYTHON_MODULE(backend)
{
    Py_Initialize();
    np::initialize();
    
    scope().attr("epsilon") = epss;
    
    class_<Point>("Point", init<>())
        .def("__len__", &Point::size)
        .def("__getitem__", static_cast<coordinate_t (Point::*)(const dimensions_t) const>(&Point::operator[]))
        .def("__str__", &Point::str)
    ;
    
    class_<Curve>("Curve", init<np::ndarray>())
        .add_property("dimensions", &Curve::dimensions)
        .add_property("points", &Curve::size)
        .def("__getitem__", static_cast<Point (Curve::*)(const curve_size_t) const>(&Curve::operator[]))
        .def("__len__", &Curve::size)
        .def("__str__", &Curve::str)
    ;
    
    class_<Curves>("Curves", init<>())
        .def("add", &Curves::add)
        .def("__getitem__", &Curves::get)
        .def("__len__", &Curves::number)
        .add_property("m", &Curves::get_m)
        .def("__str__", &Curves::str)
    ;
    
    class_<Clustering::Clustering_Result>("Clustering_Result", init<>())
        .def("__getitem__", &Clustering::Clustering_Result::get)
        .def("__len__", &Clustering::Clustering_Result::size)
        .add_property("value", &Clustering::Clustering_Result::value)
        .add_property("time", &Clustering::Clustering_Result::running_time)
        .add_property("assignment", &Clustering::Clustering_Result::assignment)
    ;
    
    class_<Clustering::Cluster_Assignment>("Clustering_Assignment", init<>())
        .def("__len__", &Clustering::Cluster_Assignment::size)
        .def("count", &Clustering::Cluster_Assignment::count)
        .def("get", &Clustering::Cluster_Assignment::get)
    ;
    
    class_<fc::Result>("Continuous_Frechet_Result", init<>())
        .add_property("time_searches", &Frechet::Continuous::Result::time_searches)
        .add_property("time_bounds", &Frechet::Continuous::Result::time_bounds)
        .add_property("number_searches", &Frechet::Continuous::Result::number_searches)
        .add_property("value", &Frechet::Continuous::Result::value)
    ;
    
    class_<fd::Result>("Discrete_Frechet_Result", init<>())
        .add_property("time", &fd::Result::time)
        .add_property("value", &fd::Result::value)
    ;
    
    class_<Coreset::Onemedian_Coreset>("Onemedian_coreset")
        .add_property("lambd", &Coreset::Onemedian_Coreset::get_lambda)
        .add_property("Lambd", &Coreset::Onemedian_Coreset::get_Lambda)
        .def("curves", &Coreset::Onemedian_Coreset::get_curves)
        .add_property("cost", &Coreset::Onemedian_Coreset::get_cost)
    ;
    
    def("dimension_reduction", jl_transform, jl_transform_overloads());
    def("continuous_frechet", continuous_frechet, continuous_frechet_overloads());
    def("discrete_frechet", discrete_frechet);
    def("discrete_kcenter", kcenter, kcenter_overloads());
    def("discrete_kmedian", kmedian, kmedian_overloads());
    def("discrete_onemedian_sampling", onemedian_sampling, onemedian_sampling_overloads());
    def("discrete_onemedian_exhaustive", onemedian_exhaustive, onemedian_exhaustive_overloads());
    def("onemedian_coreset", onemedian_coreset, onemedian_coreset_overloads());
}
