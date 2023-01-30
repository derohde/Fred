/*
Copyright 2020-2021 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>

#include <pybind11/pybind11.h>

#include "config.hpp"
#include "types.hpp"
#include "curve.hpp"
#include "frechet.hpp"

namespace py = pybind11;

namespace Simplification {

class Subcurve_Shortcut_Graph {

    Curve& curve;
    std::vector<Distances> edges;
    
public:
    
    Subcurve_Shortcut_Graph(const Curve&);
    
    Curve minimum_error_simplification(const curve_size_t) const;
};

Curve approximate_minimum_link_simplification(const Curve&, const distance_t);
Curve approximate_minimum_error_simplification(const Curve&, const curve_size_t);
 
};
