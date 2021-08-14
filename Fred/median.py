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
