# Fred ![alt text](https://raw.githubusercontent.com/derohde/Fred/master/logo/logo.png "Fred logo")
A fast, scalable and light-weight C++ Fréchet distance library, exposed to python.

## Ingredients
- continous Fréchet distance
- discrete Fréchet distance
- k-center clustering (longer running-time)
- k-median clustering (longer running-time)
- one-median clustering (short running-time)
- dimension reduction

## Installation
Get requirements under Ubuntu: `make pre`

Python3 installation into userdir: `make python3`

Python2 installation into userdir: `make python2`

## Test
Just run `python py/test.py`. Output should be: 
```
0.125
0.5
0.25
1.0
```
