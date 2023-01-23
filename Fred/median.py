"""
Copyright 2020 Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
from . import backend
from .hd_fs import curve_within_radii

import numpy as np
import itertools
import progressbar

def frechet_median(T, epsilon):
    coreset = backend.onemedian_coreset(T, epsilon/2.)
    cost = coreset.cost
    lower = cost/6
    lambd = coreset.lambd
    Lambd = coreset.Lambd
    curves = list(coreset.curves())
    curves_dist = list(set(curves))
    coreset_dist = backend.Curves()
    coreset_dist_multiplicities = [curves.count(elem) for elem in curves_dist]
    for curve_ind in curves_dist:
        coreset_dist.add(T[int(curve_ind)])
    radiis = list()
    length = 1
    for i in range(0, len(curves_dist)):
        range_ = np.linspace(0, cost, num=np.ceil(cost/(epsilon*cost/24*lambd[curves_dist[i]]/Lambd)), endpoint=True)
        length *= len(range_)
        radiis.append(range_)
    candidates = list()
    count = 0
    with progressbar.ProgressBar(max_value=length) as bar:
        for radii in itertools.product(*radiis):
            count += 1
            bar.update(count)
            rsum = np.sum(radii)
            if rsum < lower or rsum > cost:
                continue
            candidate = curve_within_radii(coreset_dist, radii)
            if candidate[0]:
                candidates.append(candidate[1])
    min_cost_curve = 0
    min_cost = np.infty
    for i in progressbar.progressbar(range(0, len(candidates))):
        cost = 0
        for j in range(0, len(coreset_dist)):
            cost += coreset_dist_multiplicities[j] * Lambd/lambd[curves_dist[j]] * 1/len(curves) * backend.continuous_frechet(candidates[i], coreset_dist[j]).value
        if cost < min_cost:
            min_cost_curve = i
            min_cost = cost
    return candidates[min_cost_curve]
