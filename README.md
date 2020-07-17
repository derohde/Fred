# Fred ![alt text](https://raw.githubusercontent.com/derohde/Fred/master/logo/logo.png "Fred logo")
A fast, scalable and light-weight C++ Fréchet distance library, exposed to python.

## Ingredients C++ Backend
`import Fred.backend as fred`

### Number of Threads

By default, Fred will automatically determine the number of threads to use. If you want to set an upper limit, call `fred.set_maximum_number_threads(number)`.

### Curve
- signature: `fred.Curve(np.ndarray)`

### Curves
- signature: `fred.Curves()`
- methods: fred.Curves.add: add curve, fred.Curves[i]: get ith curve, len(fred.Curves): number curves

#### continous Fréchet distance
- signature: `fred.continuous_frechet(curve1, curve2)`
- returns: `fred.Continuous_Frechet_Result` with members `value`, `time_bounds`: running-time for upper and lower bound, `number_searches`: number of free space diagrams built, `time_searches`: running-time for free spaces

###### continuous Frechet distance config
- approximation error: `fred.set_continuous_frechet_epsilon(epsilon)` with parameter `epsilon`, which defaults to 0.001
- rounding (rounding up to 3 decimals): `fred.set_continuous_frechet_rounding(round)` with parameter `round`, which defaults to true

#### discrete Fréchet distance
- signature: `fred.discrete_frechet(curve1, curve2)`
- returns: `fred.Discrete_Frechet_Result` with members `value` and `time`

### Clustering

#### discrete (k,l)-center clustering (continuous Fréchet)
- from [**Approximating (k,l)-center clustering for curves**](https://dl.acm.org/doi/10.5555/3310435.3310616)
- signature: `fred.discrete_klcenter(k, l, curves, with_assignment)`, `with_assignment`: defaults to false; assigns curves to nearest centers if true
- returns: `fred.Clustering_Result` with mebers `value`: objective value, `time`: running-time, `assignment`: empty if with_assignment=false

#### discrete (k,l)-median clustering (continuous Fréchet)
- Algorithm 6 in [**Coresets for (k,l)-Clustering under the Fréchet distance**](https://arxiv.org/pdf/1901.01870.pdf) + simplification
- signature: `fred.discrete_klmedian(k, l, curves, with_assignment)` with parameter `with_assignment`: defaults to false; assigns curves to nearest centers if true
- returns: `fred.Clustering_Result` with mebers `value`: objective value, `time`: running-time, `assignment`: empty if with_assignment=false

#### discrete one-median clustering (continuous Fréchet) via sampling 
- [Section 3 in **Random Projections and Sampling Algorithms for Clustering of High Dimensional Polygonal Curves**](https://papers.nips.cc/paper/9443-random-projections-and-sampling-algorithms-for-clustering-of-high-dimensional-polygonal-curves)
- signature: `fred.discrete_onemedian_sampling(curves, epsilon_sampling, with_assignment)` with parameters `epsilon_sampling`: (1+epsilon) approximation parameter, `with_assignment`: defaults to false; assigns curves to nearest centers if true
- returns: `fred.Clustering_Result` with mebers `value`: objective value, `time`: running-time, `assignment`: empty if with_assignment=false

#### Clustering Result
- signature: `fred.Clustering_Result`
- methods: `len(fred.Clustering_Result)`: number of centers, `fred.Clustering_Result[i]`: get ith center
- members: `value`: objective value, `time`: running-time, `assignment`: empty if with_assignment=false

#### Cluster Assignment
- signature: `fred.Cluster_Assignment`
- methods: len(fred.Cluster_Assignment): number of centers, fred.Cluster_Assignment.count(i): number of curves assigned to center i, fred.Cluster_Assignment.get(i,j): get index of jth curve assigned to center i

### dimension reduction via. gaussian random projection 
- [Section 2 in **Random Projections and Sampling Algorithms for Clustering of High Dimensional Polygonal Curves**](https://papers.nips.cc/paper/9443-random-projections-and-sampling-algorithms-for-clustering-of-high-dimensional-polygonal-curves)
- signature: `fred.dimension_reduction(curves, epsilon, empirical_constant)` with parameters `epsilon`: (1+epsilon) approximation parameter, `empirical_constant`: use constant of empirical study (faster, but less accurate)
- returns: `fred.Curves` collection of curves
  
## Mini Example
```python
import Fred.backend as fred
import numpy as np
import pandas as pd

curve1d = fred.Curve(np.array([1., 2.])) # Curve stores a polygonal curve with 
                                         # at least two points of at least one 
                                         # and equal number of dimensions

curve2d1 = fred.Curve(np.array([[1., 0.], [2., 0.], [3., 0.]])) # any number of dimensions and points works
curve2d2 = fred.Curve(np.array([[1., -1.], [2., -1.], [3., -1.]])) 

print(curve2d1)

print("distance is {}".format(fred.continuous_frechet(curve2d1, curve2d2).value))

print("download HUGE curves")

import requests, zipfile, io             # you can use all libraries 
                                         # that work with numpy to read data into fred
                                         
re = requests.get("https://archive.ics.uci.edu/ml/machine-learning-databases/00447/data.zip", stream=True)
zf = zipfile.ZipFile(io.BytesIO(re.content))

ps1 = fred.Curve(pd.read_csv(zf.open('PS1.txt'), delimiter="\t", header=None).values)
ps2 = fred.Curve(pd.read_csv(zf.open('PS2.txt'), delimiter="\t", header=None).values)
ps3 = fred.Curve(pd.read_csv(zf.open('PS3.txt'), delimiter="\t", header=None).values)
ps4 = fred.Curve(pd.read_csv(zf.open('PS4.txt'), delimiter="\t", header=None).values)
ps5 = fred.Curve(pd.read_csv(zf.open('PS5.txt'), delimiter="\t", header=None).values)
ps6 = fred.Curve(pd.read_csv(zf.open('PS6.txt'), delimiter="\t", header=None).values)

curves = fred.Curves() # for clustering or if you want to apply dimension reduction
                       # you need to encapsulate your curves in a Curves object
              
curves.add(ps1)
curves.add(ps2)
curves.add(ps3)
curves.add(ps4)
curves.add(ps5)
curves.add(ps6)

curves = fred.dimension_reduction(curves, 0.95) # fred is pretty fast but with high dimensional data
                                                # a dimension reduction massively improves running-time
                                                # even for smaller values of epsilon
                                  
                                               
clustering = fred.discrete_klmedian(2,, 100, curves)

print("clustering cost is {}".format(clustering.value))

for i, center in enumerate(clustering):
    print("center {} is {}".format(i, center))
```
  
## Installation
Get requirements under Ubuntu: `make pre`

Python3 installation into userdir: `make python3`

Python2 installation into userdir: `make python2`

## Test
Just run `python py/test.py`.
