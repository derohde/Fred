import Fred.backend as fred
import Fred
import numpy as np


nzigzags = 20
hsep = 10
vsep = 20


curves = fred.Curves()

p1 = np.empty((nzigzags, 2))
p2 = np.empty((nzigzags, 2))

for i in range(nzigzags):
    p1[i, 0] = 0 + i * hsep/(nzigzags-1)
    p2[i, 1] = 0 + i * vsep/(nzigzags-1)
    if i % 2 == 0:
        p1[i, 1] = 0
        p2[i, 0] = 0
    else:
        p1[i, 1] = vsep
        p2[i, 0] = hsep

curves.add(fred.Curve(p1, "Input Curve {}".format(1)))
curves.add(fred.Curve(p2, "Input Curve {}".format(2)))

clustering = fred.two_two_dtw_one_two_median(curves)
clustering_e = fred.two_two_dtw_one_two_median_exact(curves)

Fred.plot_curve(curves, clustering, clustering_e)
