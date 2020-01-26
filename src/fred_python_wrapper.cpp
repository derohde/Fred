/*
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <boost/python.hpp>

#include "frechet.hpp"
#include "curve.hpp"
#include "jl_transform.hpp"
#include "clustering.hpp"

using namespace boost::python;
namespace np = boost::python::numpy;
namespace fc = Frechet::Continuous;
namespace fd = Frechet::Discrete;

distance_t epss = 0.001;

fc::Result continuous_frechet(const Curve &curve1, const Curve &curve2, const distance_t eps) { 
    return fc::distance(curve1, curve2, eps);
}

fd::Result discrete_frechet(const Curve &curve1, const Curve &curve2) {
    return fd::distance(curve1, curve2);
}

Curves jl_transform(const Curves &in, const double epsilon, const bool empirical_constant) {
    
    Curves curvesrp = JLTransform::transform_naive(in, epsilon, empirical_constant);
    
    return curvesrp;
}

Clustering::Clustering_Result kcenter(const curve_size_t num_centers, const Curves &in, const distance_t eps, const bool with_assignment = false) {
    auto result = Clustering::gonzalez(num_centers, in, eps, false, with_assignment);
    
    return result;
}

Clustering::Clustering_Result onemedian_sampling(const Curves &in, const double epsilon, const distance_t eps, const bool with_assignment = false) {
    
    auto result = Clustering::one_median_sampling(epsilon, in, eps, with_assignment);
    
    return result;
}

Clustering::Clustering_Result onemedian_exhaust(const Curves &in, const distance_t eps, const bool with_assignment = false) {

    auto result = Clustering::one_median_exhaustive(in, eps, with_assignment);
    
    return result;
}

Clustering::Clustering_Result kmedian(const curve_size_t num_centers, const Curves &in, const distance_t eps, const bool with_assignment = false) {

    auto result = Clustering::arya(num_centers, in, eps, with_assignment);
    
    return result;
}

BOOST_PYTHON_MODULE(Fred)
{
    Py_Initialize();
    np::initialize();
    
    scope().attr("epsilon") = epss;
    
    class_<Curve>("Curve", init<np::ndarray>())
        .add_property("dimensions", &Curve::dimensions)
        .add_property("points", &Curve::size)
        .def("__len__", &Curve::size)
    ;

    class_<Curves>("Curves", init<>())
        .def("add", static_cast<void (Curves::*)(const Curve&)>(&Curves::push_back))
        .def("__getitem__", &Curves::get)
        .def("__len__", &Curves::size)
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
        .def("__getitem__", &Clustering::Cluster_Assignment::get)
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
    
    
    def("dimension_reduction", jl_transform);
    def("continuous_frechet", continuous_frechet);
    def("discrete_frechet", discrete_frechet);
    def("discrete_kcenter", kcenter);
    def("discrete_kmedian", kmedian);
    def("discrete_onemedian_sampling", onemedian_sampling);
    def("discrete_onemedian_exhaustive", onemedian_exhaust);
}
