# Fred ![alt text](https://raw.githubusercontent.com/derohde/Fred/master/logo/logo.png "Fred logo")
A fast, scalable and light-weight C++ Fréchet distance library, exposed to python and focused on (k,l)-clustering of polygonal curves.

### NOW USING PYBIND11 INSTEAD OF BOOST!
### NOW AVAILABLE VIA PIP

## Ingredients C++ Backend
`import Fred.backend as fred`

### Number of Threads

By default, Fred will automatically determine the number of threads to use. If you want to set an upper limit, call `fred.set_maximum_number_threads(number)`.

### Curve
- signature: `fred.Curve(np.ndarray)`, `fred.Curve(np.ndarray, str name)`
- properties: `fred.Curve.values`: curves as `np.ndarray`, `fred.Curve.name`: get name of curve, `fred.Curve.dimensions`: dimension of curve, `fred.Curve.complexity`: number of points of curve

### Curves
- signature: `fred.Curves()`
- methods: `fred.Curves.add(curve)`: add curve, `fred.Curves[i]`: get ith curve, `len(fred.Curves)`: number curves, `fred.Curves.simplify(l)`: return set of simplified curves
- properties:  `fred.Curves.m`: maximum complexity of the contained curves, `fred.Curves.values`: curves as `np.ndarray`

#### continous Fréchet distance
- signature: `fred.continuous_frechet(curve1, curve2)`
- returns: `fred.Continuous_Frechet_Result` with members `value`, `time_bounds`: running-time for upper and lower bound, `number_searches`: number of free space diagrams built, `time_searches`: running-time for free spaces

###### continuous Frechet distance config
- approximation error: `fred.set_continuous_frechet_epsilon(epsilon)` with parameter `epsilon`, which defaults to 0.001
- rounding (rounding up to 3 decimals): `fred.set_continuous_frechet_rounding(round)` with parameter `round`, which defaults to true

#### discrete Fréchet distance
- signature: `fred.discrete_frechet(curve1, curve2)`
- returns: `fred.Discrete_Frechet_Result` with members `value` and `time`

#### discrete dynamic time warping distance
- signature: `fred.discrete_dynamic_time_warping(curve1, curve2)`
- returns: `fred.Discrete_Dynamic_Time_Warping_Distance` with members `value` and `time`

### Curve Simplification

#### weak minimum error simplification
- graph approach from [**Polygonal Approximations of a Curve — Formulations and Algorithms**](https://www.sciencedirect.com/science/article/pii/B9780444704672500114)
- signature: `fred.weak_minimum_error_simplification(fred.Curve, int complexity)`
- returns: `fred.Curve`that uses input curves vertices, with `complexity` number of vertices and that has minimum distance to input curve

#### approximate weak minimum link simplification
- algorithm "FS" from [**Near-Linear Time Approximation Algorithms for Curve Simplification**](https://link.springer.com/article/10.1007/s00453-005-1165-y)
- signature: `fred.approximate_weak_minimum_error_simplification(fred.Curve, double error)`
- returns: `fred.Curve` that uses input curves vertices, is of small complexity and with distance to input curve at most `error`

#### approximate weak minimum error simplification
- binary search on `fred.approximate_weak_minimum_link_simplification`
- signature: `fred.approximate_weak_minimum_error_simplification(fred.Curve, int complexity)`
- returns: `fred.Curve`that uses input curves vertices, with `complexity` number of vertices and that has small distance to input curve

### Clustering

##### Distance_Matrix

A `fred.Distance_Matrix()` can be used to speed up consecutive calls of `fred.discrete_klcenter` and `fred.discrete_klmedian`. As the name suggests, it stores the distances already computed.

#### discrete (k,l)-center clustering (continuous Fréchet)
- from [**Approximating (k,l)-center clustering for curves**](https://dl.acm.org/doi/10.5555/3310435.3310616)
- signature: `fred.discrete_klcenter_multi(k, l, curves, distances, center_domain, random_first_center)` with parameters 
    - `k`: number of centers
    - `l`: maximum complexity of the centers, only used when center_domain is default value
    - `distances`: `fred.Distance_Matrix`, defaults to empty `fred.Distance_Matrix`
    - `center_domain`: possible centers, defaults to empty `fred.Curves()`, in this case the input is simplified and used as center domain
    - `random_first_center`: determines if first center is chosen uniformly at random or first curve is used as first center, optional, defaults to true
- returns: `fred.Clustering_Result` with mebers 
    - `value`: objective value 
    - `time`: running-time 
    - `assignment`: empty if compute_assignment has not been called

#### discrete (k,l)-median clustering (continuous Fréchet)
- Algorithm from section 4.3 in [**Geometric Approximation Algorithms**](http://www.ams.org/books/surv/173/) + simplification
- signature: `fred.discrete_klmedian_multi(k, l, curves, distances, center_domain)` with parameters 
    - `k`: number of centers
    - `l`: maximum complexity of the centers, only used when center_domain is default value
    - `distances`: `fred.Distance_Matrix`, defaults to empty `fred.Distance_Matrix`
    - `center_domain`: possible centers, optional parameter, if not given the input is simplified and used as center domain
- returns: `fred.Clustering_Result` with mebers 
    - `value`: objective value 
    - `time`: running-time 
    - `assignment`: empty if compute_assignment has not been called

#### Clustering Result
- signature: `fred.Clustering_Result`
- methods: `len(fred.Clustering_Result)`: number of centers, `fred.Clustering_Result[i]`: get ith center, `fred.Clustering_Result.compute_assignment(fred.Curves)`: assigns every curve to its nearest center
- members: `value`: objective value, `time`: running-time, `assignment`: empty if compute_assignment was not called

#### Cluster Assignment
- signature: `fred.Cluster_Assignment`
- methods: `len(fred.Cluster_Assignment)`: number of centers, `fred.Cluster_Assignment.count(i)`: number of curves assigned to center i, `fred.Cluster_Assignment.get(i,j)`: get index of jth curve assigned to center i

### Dimension Reduction via Gaussian Random Projection 
- [Section 2 in **Random Projections and Sampling Algorithms for Clustering of High Dimensional Polygonal Curves**](https://papers.nips.cc/paper/9443-random-projections-and-sampling-algorithms-for-clustering-of-high-dimensional-polygonal-curves)
- signature: `fred.dimension_reduction(curves, epsilon, empirical_constant)` with parameters `epsilon`: (1+epsilon) approximation parameter, `empirical_constant`: use constant of empirical study (faster, but less accurate), defaults to `True`
- returns: `fred.Curves` collection of curves
  
## Installation

### Requirements

You have to have installed:
 - git
 - openmp available (should be a part of your compiler)
 
Thats it!

### Installation Procedure

 - Variant 1: simply run `pip install Fred-Frechet`
 - Variant 2: clone repository and run `make` for installation into userdir

## Test
Just run `python py/test.py`.
  
## Mini Example
```python
import Fred.backend as fred
import Fred
import numpy as np
import pandas as pd

curve1d = fred.Curve(np.array([1., 2.])) # Curve stores a polygonal curve with 
                                         # at least two points of at least one 
                                         # and equal number of dimensions

curve2d1 = fred.Curve(np.array([[1., 0.], [2., 1.], [3., 0.]])) # any number of dimensions and points works
curve2d2 = fred.Curve(np.array([[1., -1.], [2., -2.], [3., -1.]]), "optional name, e.g. displayed in plot") 

print(curve2d1)

Fred.plot_curve(curve2d1, curve2d2)
Fred.plot_curve(curve2d2, fred.weak_minimum_error_simplification(curve2d2, 2))

print("distance is {}".format(fred.continuous_frechet(curve2d1, curve2d2).value))

print("download HUGE curves") 

import requests, zipfile, io             # you can use all libraries 
                                         # that work with numpy to read data into fred
                                         
re = requests.get("https://archive.ics.uci.edu/ml/machine-learning-databases/00447/data.zip", stream=True)
zf = zipfile.ZipFile(io.BytesIO(re.content))

ps1 = fred.Curve(pd.read_csv(zf.open('PS1.txt'), delimiter="\t", header=None).values[:50], "PS1")
ps2 = fred.Curve(pd.read_csv(zf.open('PS2.txt'), delimiter="\t", header=None).values[:50], "PS2")
ps3 = fred.Curve(pd.read_csv(zf.open('PS3.txt'), delimiter="\t", header=None).values[:50], "PS3")
ps4 = fred.Curve(pd.read_csv(zf.open('PS4.txt'), delimiter="\t", header=None).values[:50], "PS4")
ps5 = fred.Curve(pd.read_csv(zf.open('PS5.txt'), delimiter="\t", header=None).values[:50], "PS5")
ps6 = fred.Curve(pd.read_csv(zf.open('PS6.txt'), delimiter="\t", header=None).values[:50], "PS6")

curves = fred.Curves() # for clustering or if you want to apply dimension reduction
                       # you need to encapsulate your curves in a Curves object
              
curves.add(ps1)
curves.add(ps2)
curves.add(ps3)
curves.add(ps4)
curves.add(ps5)
curves.add(ps6)

Fred.plot_curve(curves)

curves = fred.dimension_reduction(curves, 0.95) # fred is pretty fast but with high dimensional data
                                                # a dimension reduction massively improves running-time
                                                # even for smaller values of epsilon
                                                
Fred.plot_curve(curves)
                                  
# Oneshot clustering - if you already know the value of k
                                  
clustering = fred.discrete_klcenter(2, 10, curves) # fast but coarse
          
clustering = fred.discrete_klmedian(2, 10, curves) # slow but better results

print("clustering cost is {}".format(clustering.value))

for i, center in enumerate(clustering):
    print("center {} is {}".format(i, center))
    
    
Fred.plot_curve(clustering)

# Multiple clustering calls - if you need to find a suitable value for k

dm = fred.Distance_Matrix() # computing the Fréchet distance is costly,
                            # therefore we buffer each distance already computed to
                            # speed up consecutive clustering calls
                            
for k in range(2, 6):
    
    clustering = fred.discrete_klcenter(k, 10, curves, dm)
    print("clustering cost is {}".format(clustering.value))
            
    clustering = fred.discrete_klmedian(k, 10, curves, dm)
    print("clustering cost is {}".format(clustering.value))
    
clustering.compute_assignment(curves)

for i in range(0, len(clustering)):
    for j in range(0, clustering.assignment.count(i)):
        print("{} was assigned to center {}".format(curves[clustering.assignment.get(i,j)].name, clustering[i].name))
```
