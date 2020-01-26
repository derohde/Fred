# Fred ![alt text](https://raw.githubusercontent.com/derohde/Fred/master/logo/logo.png "Fred logo")
A fast, scalable and light-weight C++ Fréchet distance library, exposed to python.

## Ingredients
- continous Fréchet distance
  - signature: `Fred.continuous_frechet(curve1, curve2, approximation_error)` default for approximation_error is 0.001
  - returns: `Fred.Continuous_Frechet_Result` with members `value`, `time_bounds`: running-time for upper and lower bound, `number_searches`: number of free space diagrams built, `time_searches`: running-time for free spaces
- discrete Fréchet distance
  - signature: `Fred.discrete_frechet(curve1, curve2)`
  - returns: `Fred.Discrete_Frechet_Result` with members `value` and `time`
- discrete k-center clustering (continuous Fréchet) [Without simplification; from **Approximating (k,l)-center clustering for curves**](https://dl.acm.org/doi/10.5555/3310435.3310616)
  - signature: `Fred.discrete_kcenter(k, curves, approximation_error, with_assignment)` with parameters `approximation_error`: see continuous Fréchet, `with_assignment`: defaults to false; assigns curves to nearest centers if true
  - returns: `Fred.Clustering_Result` with mebers `value`: objective value, `time`, `assignment`: empty if with_assignment=false
- discrete k-median clustering (continuous Fréchet) [Algorithm 6 in **Coresets for (k,l)-Clustering under the Fréchet distance**](https://arxiv.org/pdf/1901.01870.pdf)
  - signature: `Fred.discrete_kmedian(k, curves, approximation_error, with_assignment)` with parameters `approximation_error`: see continuous Fréchet, `with_assignment`: defaults to false; assigns curves to nearest centers if true
  - returns: `Fred.Clustering_Result` with mebers `value`: objective value, `time`, `assignment`: empty if with_assignment=false
- discrete one-median clustering (continuous Fréchet) via sampling [Section 3 in **Random Projections and Sampling Algorithms for Clustering of High Dimensional Polygonal Curves**](https://papers.nips.cc/paper/9443-random-projections-and-sampling-algorithms-for-clustering-of-high-dimensional-polygonal-curves)
  - signature: `Fred.discrete_onemedian_sampling(curves, epsilon_sampling, approximation_error, with_assignment)` with parameters `approximation_error`: see continuous Fréchet, `epsilon_sampling`: (1+epsilon) approximation parameter, `with_assignment`: defaults to false; assigns curves to nearest centers if true
  - returns: `Fred.Clustering_Result` with mebers `value`: objective value, `time`, `assignment`: empty if with_assignment=false
- dimension reduction via. gaussian random projection [Section 2 in **Random Projections and Sampling Algorithms for Clustering of High Dimensional Polygonal Curves**](https://papers.nips.cc/paper/9443-random-projections-and-sampling-algorithms-for-clustering-of-high-dimensional-polygonal-curves)
  - signature: `Fred.dimension_reduction(curves, epsilon, empirical_constant)` with parameters `epsilon`: (1+epsilon) approximation parameter, `empirical_constant`: use constant of empirical study (faster, but less accurate)
  - returns: `Fred.Curves` collection of curves
  
## Installation
Get requirements under Ubuntu: `make pre`

Python3 installation into userdir: `make python3`

Python2 installation into userdir: `make python2`

## Test
Just run `python py/test.py`.
