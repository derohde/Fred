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
#include "simplification.hpp"
#include "dynamic_time_warping.hpp"

using namespace boost::python;
namespace np = boost::python::numpy;
namespace fc = Frechet::Continuous;
namespace fd = Frechet::Discrete;
namespace ddtw = Dynamic_Time_Warping::Discrete;

const distance_t default_epsilon = 0.001;

fc::Distance continuous_frechet(const Curve &curve1, const Curve &curve2) { 
    return fc::distance(curve1, curve2);
}

fd::Distance discrete_frechet(const Curve &curve1, const Curve &curve2) {
    return fd::distance(curve1, curve2);
}

ddtw::Distance discrete_dynamic_time_warping(const Curve &curve1, const Curve &curve2) {
    return ddtw::distance(curve1, curve2);
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

Clustering::Clustering_Result dtw_one_median(const Curves &in) {
    auto result = Clustering::dtw_one_median(in);
    return result;
}

Clustering::Clustering_Result klcenter_multi(const curve_number_t num_centers, const curve_size_t ell, const Curves &in, Clustering::Distance_Matrix &distances,
                                       bool with_assignment = false, const Curves &center_domain = Curves()) {
    auto result = Clustering::gonzalez(num_centers, ell, in, distances, false, with_assignment, center_domain);
    return result;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(klcenter_multi_overloads, klcenter_multi, 4, 6);

Clustering::Clustering_Result klcenter(const curve_number_t num_centers, const curve_size_t ell, const Curves &in, bool with_assignment = false, const Curves &center_domain = Curves()) {
    Clustering::Distance_Matrix distances;
    auto result = Clustering::gonzalez(num_centers, ell, in, distances, false, with_assignment, center_domain);
    return result;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(klcenter_overloads, klcenter, 3, 5);

Clustering::Clustering_Result onemedian_sampling(const curve_size_t ell,  Curves &in, const double epsilon, const bool with_assignment = false, const Curves &center_domain = Curves()) {
    
    auto result = Clustering::one_median_sampling(ell, in, epsilon, with_assignment);
    
    return result;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(onemedian_sampling_overloads, onemedian_sampling, 3, 5);

Clustering::Clustering_Result onemedian_exhaustive(const curve_size_t ell,  Curves &in, const bool with_assignment = false, const Curves &center_domain = Curves()) {

    auto result = Clustering::one_median_exhaustive(ell, in, with_assignment);
    
    return result;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(onemedian_exhaustive_overloads, onemedian_exhaustive, 2, 4);

Clustering::Clustering_Result klmedian_multi(const curve_number_t num_centers, const curve_size_t ell, const Curves &in, Clustering::Distance_Matrix distances,
                                             bool with_assignment = false, const Curves &center_domain = Curves()) {

    auto result = Clustering::arya(num_centers, ell, in, distances, with_assignment, center_domain);
    
    return result;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(klmedian_multi_overloads, klmedian_multi, 4, 6);

Clustering::Clustering_Result klmedian(const curve_number_t num_centers, const curve_size_t ell, const Curves &in, const bool with_assignment = false, const Curves &center_domain = Curves()) {

    auto distances = Clustering::Distance_Matrix();
    auto result = Clustering::arya(num_centers, ell, in, distances, with_assignment, center_domain);
    
    return result;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(klmedian_overloads, klmedian, 3, 5);

Coreset::Onemedian_Coreset onemedian_coreset(const Curves &in, const curve_size_t ell, const double epsilon, const double constant = 1) {
    return Coreset::Onemedian_Coreset(ell, in, epsilon, constant);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(onemedian_coreset_overloads, onemedian_coreset, 3, 4);

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

BOOST_PYTHON_MODULE(backend)
{
    Py_Initialize();
    np::initialize();
    
    scope().attr("default_epsilon_continuous_frechet") = default_epsilon;
    
    class_<Point>("Point", init<dimensions_t>())
        .def("__len__", &Point::dimensions)
        .def("__getitem__", &Point::get)
        .def("__str__", &Point::str)
        .def("__iter__", iterator<Point>())
        .def("__repr__", &Point::repr)
        .add_property("values", &Point::as_ndarray)
    ;
    
    class_<Points>("Points", init<dimensions_t>())
        .def("__len__", &Points::number)
        .def("__getitem__", &Points::get)
        .def("__str__", &Points::str)
        .def("__iter__", iterator<Points>())
        .def("__repr__", &Points::repr)
        .add_property("values", &Points::as_ndarray)
        .add_property("centroid", &Points::centroid)
    ;
    
    class_<Curve>("Curve", init<np::ndarray>())
        .def(init<np::ndarray, std::string>())
        .add_property("dimensions", &Curve::dimensions)
        .add_property("complexity", &Curve::complexity)
        .add_property("name", &Curve::get_name)
        .add_property("values", &Curve::as_ndarray)
        .add_property("centroid", &Curve::centroid)
        .def("__getitem__", &Curve::get)
        .def("__len__", &Curve::complexity)
        .def("__str__", &Curve::str)
        .def("__iter__", iterator<Curve>())
        .def("__repr__", &Curve::repr)
    ;
        
    class_<Curves>("Curves", init<>())
        .add_property("m", &Curves::get_m)
        .def("add", &Curves::add)
        .def("simplify", &Curves::simplify)
        .def("__getitem__", &Curves::get)
        .def("__len__", &Curves::number)
        .def("__str__", &Curves::str)
        .def("__iter__", iterator<Curves>())
        .def("__repr__", &Curves::repr)
    ;
    
    class_<fc::Distance>("Continuous_Frechet_Distance", init<>())
        .add_property("time_searches", &fc::Distance::time_searches)
        .add_property("time_bounds", &fc::Distance::time_bounds)
        .add_property("number_searches", &fc::Distance::number_searches)
        .add_property("value", &fc::Distance::value)
        .def("__repr__", &fc::Distance::repr)
    ;
    
    class_<fd::Distance>("Discrete_Frechet_Distance", init<>())
        .add_property("time", &fd::Distance::time)
        .add_property("value", &fd::Distance::value)
        .def("__repr__", &fd::Distance::repr)
    ;
    
    class_<ddtw::Distance>("Discrete_Dynamic_Time_Warping_Distance", init<>())
        .add_property("time", &ddtw::Distance::time)
        .add_property("value", &ddtw::Distance::value)
        .def("__repr__", &ddtw::Distance::repr)
    ;
    
    class_<Clustering::Distance_Matrix>("Distance_Matrix", init<>());
    
    class_<Clustering::Clustering_Result>("Clustering_Result", init<>())
        .add_property("value", &Clustering::Clustering_Result::value)
        .add_property("time", &Clustering::Clustering_Result::running_time)
        .add_property("assignment", &Clustering::Clustering_Result::assignment)
        .def("__getitem__", &Clustering::Clustering_Result::get)
        .def("__len__", &Clustering::Clustering_Result::size)
        .def("__iter__", range(&Clustering::Clustering_Result::cbegin, &Clustering::Clustering_Result::cend))
    ;
    
    class_<Clustering::Cluster_Assignment>("Cluster_Assignment", init<>())
        .def("__len__", &Clustering::Cluster_Assignment::size)
        .def("count", &Clustering::Cluster_Assignment::count)
        .def("get", &Clustering::Cluster_Assignment::get)
    ;
    
    /*class_<Coreset::Onemedian_Coreset>("Onemedian_coreset")
        .add_property("lambd", &Coreset::Onemedian_Coreset::get_lambda)
        .add_property("Lambd", &Coreset::Onemedian_Coreset::get_Lambda)
        .add_property("cost", &Coreset::Onemedian_Coreset::get_cost)
        .def("curves", &Coreset::Onemedian_Coreset::get_curves)
    ;*/
    
    def("set_continuous_frechet_epsilon", set_frechet_epsilon);
    def("set_continuous_frechet_rounding", set_frechet_rounding);
    def("get_continuous_frechet_epsilon", get_frechet_epsilon);
    def("get_continuous_frechet_rounding", get_frechet_rounding);
    
    def("continuous_frechet", continuous_frechet);
    def("discrete_frechet", discrete_frechet);
    def("discrete_dynamic_time_warping", discrete_dynamic_time_warping);
    
    def("weak_minimum_error_simplification", weak_minimum_error_simplification);
    def("approximate_weak_minimum_link_simplification", approximate_weak_minimum_link_simplification);
    def("approximate_weak_minimum_error_simplification", approximate_weak_minimum_error_simplification);
    
    def("dimension_reduction", jl_transform, jl_transform_overloads());

    def("discrete_klcenter", klcenter, klcenter_overloads());
    def("discrete_klmedian", klmedian, klmedian_overloads());
    def("discrete_klcenter_multi", klcenter_multi, klcenter_multi_overloads());
    def("discrete_klmedian_multi", klmedian_multi, klmedian_multi_overloads());
    
    // these are experimental
    def("dtw_one_median", dtw_one_median);
    def("discrete_onemedian_sampling", onemedian_sampling, onemedian_sampling_overloads());
    def("discrete_onemedian_exhaustive", onemedian_exhaustive, onemedian_exhaustive_overloads());
    def("onemedian_coreset", onemedian_coreset, onemedian_coreset_overloads());
    
    
    def("set_maximum_number_threads", set_number_threads);
}
