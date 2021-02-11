# Fred ![alt text](https://raw.githubusercontent.com/derohde/Fred/master/logo/logo.png "Fred logo")
A fast, scalable and light-weight C++ Fréchet distance library, exposed to python and focused on (k,l)-clustering of polygonal curves.

## Ingredients C++ Backend
`import Fred.backend as fred`

### Number of Threads

By default, Fred will automatically determine the number of threads to use. If you want to set an upper limit, call `fred.set_maximum_number_threads(number)`.

### Curve
- signature: `fred.Curve(np.ndarray)`
- properties: `fred.Curve.values`: curves as `np.ndarray`, `fred.Curve.name`: get name of curve, `fred.Curve.dimensions`: dimension of curve, `fred.Curve.complexity`: number of points of curve

### Curves
- signature: `fred.Curves()`
- methods: `fred.Curves.add(curve)`: add curve, `fred.Curves[i]`: get ith curve, `len(fred.Curves)`: number curves, `fred.Curves.simplify(l)`: return set of simplified curves
- properties:  `fred.Curves.m`: maximum complexity of the contained curves

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

##### Distance_Matrix

A `fred.Distance_Matrix()` can be used to speed up consecutive calls of `fred.discrete_klcenter` and `fred.discrete_klmedian`. As the name suggests, it stores the Frechet distances already computed.


#### discrete (k,l)-center clustering (continuous Fréchet) -- multiple calls
- from [**Approximating (k,l)-center clustering for curves**](https://dl.acm.org/doi/10.5555/3310435.3310616)
- signature: `fred.discrete_klcenter_multi(k, l, curves, distances, with_assignment, center_domain)`with parameters `distances`: `fred.Distance_Matrix`, `with_assignment`: defaults to false; assigns curves to nearest centers if true, `center_domain`: possible centers, defaults to empty `fred.Curves()`, in this case the input is simplified and used as center domain
- returns: `fred.Clustering_Result` with mebers `value`: objective value, `time`: running-time, `assignment`: empty if with_assignment=false

#### discrete (k,l)-median clustering (continuous Fréchet) -- multiple calls
- Algorithm 6 in [**Coresets for (k,l)-Clustering under the Fréchet distance**](https://arxiv.org/pdf/1901.01870.pdf) + simplification
- signature: `fred.discrete_klmedian_multi(k, l, curves, distances, with_assignment, center_domain)` with parameters `distances`: `fred.Distance_Matrix`, `with_assignment`: defaults to false; assigns curves to nearest centers if true, `center_domain`: possible centers, defaults to empty `fred.Curves()`, in this case the input is simplified and used as center domain
- returns: `fred.Clustering_Result` with mebers `value`: objective value, `time`: running-time, `assignment`: empty if with_assignment=false


#### discrete (k,l)-center clustering (continuous Fréchet) -- oneshot
- from [**Approximating (k,l)-center clustering for curves**](https://dl.acm.org/doi/10.5555/3310435.3310616)
- signature: `fred.discrete_klcenter(k, l, curves, with_assignment, center_domain)` with parameters `with_assignment`: defaults to false; assigns curves to nearest centers if true, `center_domain`: possible centers, defaults to empty `fred.Curves()`, in this case the input is simplified and used as center domain
- returns: `fred.Clustering_Result` with mebers `value`: objective value, `time`: running-time, `assignment`: empty if with_assignment=false

#### discrete (k,l)-median clustering (continuous Fréchet) -- oneshot
- Algorithm 6 in [**Coresets for (k,l)-Clustering under the Fréchet distance**](https://arxiv.org/pdf/1901.01870.pdf) + simplification
- signature: `fred.discrete_klmedian(k, l, curves, with_assignment, center_domain)` with parameters `with_assignment`: defaults to false; assigns curves to nearest centers if true, `center_domain`: possible centers, defaults to empty `fred.Curves()`, in this case the input is simplified and used as center domain
- returns: `fred.Clustering_Result` with mebers `value`: objective value, `time`: running-time, `assignment`: empty if with_assignment=false

#### Clustering Result
- signature: `fred.Clustering_Result`
- methods: `len(fred.Clustering_Result)`: number of centers, `fred.Clustering_Result[i]`: get ith center
- members: `value`: objective value, `time`: running-time, `assignment`: empty if with_assignment=false

#### Cluster Assignment
- signature: `fred.Cluster_Assignment`
- methods: `len(fred.Cluster_Assignment)`: number of centers, `fred.Cluster_Assignment.count(i)`: number of curves assigned to center i, `fred.Cluster_Assignment.get(i,j)`: get index of jth curve assigned to center i

### Dimension Reduction via Gaussian Random Projection 
- [Section 2 in **Random Projections and Sampling Algorithms for Clustering of High Dimensional Polygonal Curves**](https://papers.nips.cc/paper/9443-random-projections-and-sampling-algorithms-for-clustering-of-high-dimensional-polygonal-curves)
- signature: `fred.dimension_reduction(curves, epsilon, empirical_constant)` with parameters `epsilon`: (1+epsilon) approximation parameter, `empirical_constant`: use constant of empirical study (faster, but less accurate)
- returns: `fred.Curves` collection of curves
  
## Installation
Get requirements under Ubuntu: `make pre`

Python3 installation into userdir: `make install`

### If something does not work with Boost

Manual installation of Boost

- `mkdir $HOME/boost` (This folder is hardcoded in setup.py, another location won't work.)
- `cd /tmp`
- `wget https://dl.bintray.com/boostorg/release/1.73.0/source/boost_1_73_0.tar.gz`
- `tar -xzf boost_1_73_0.tar.gz`
- `cd boost_1_73_0`
- `./bootstrap.sh --with-python=/usr/bin/python3`
- `./b2 install --prefix=$HOME/boost`

After that, go back to Freds folder and run `make clean` and then `make install`

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
curve2d2 = fred.Curve(np.array([[1., -1.], [2., -2.], [3., -1.]])) 

print(curve2d1)

Fred.plot_curve(curve2d1, curve2d2)
Fred.plot_curve(curve2d2, fred.weak_minimum_error_simplification(curve2d2, 2))

print("distance is {}".format(fred.continuous_frechet(curve2d1, curve2d2).value))

print("download HUGE curves")  # WARNING: running the algorithms with the following input may take several hours,
                               #          depending on your hardware

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

Fred.plot_curve(curves)

curves = fred.dimension_reduction(curves, 0.95) # fred is pretty fast but with high dimensional data
                                                # a dimension reduction massively improves running-time
                                                # even for smaller values of epsilon
                                                
Fred.plot_curve(curves)
                                  
# Oneshot clustering - if you already know the value of k
                                  
clustering = fred.discrete_klcenter(2, 100, curves) # fast but coarse
          
clustering = fred.discrete_klmedian(2, 100, curves) # slow but better results

print("clustering cost is {}".format(clustering.value))

for i, center in enumerate(clustering):
    print("center {} is {}".format(i, center))
    
    
Fred.plot_curve(clustering)

# Multiple clustering calls - if you need to find a suitable value for k

dm = fred.Distance_Matrix() # computing the Fréchet distance is costly,
                            # therefore we buffer each distance already computed to
                            # speed up consecutive clustering calls
                            
for k in range(2, 6):
    
    clustering = fred.discrete_klcenter_multi(k, 100, curves, dm)
    print("clustering cost is {}".format(clustering.value))
            
    clustering = fred.discrete_klmedian_multi(k, 100, curves, dm)
    print("clustering cost is {}".format(clustering.value))
```
