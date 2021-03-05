import Fred.backend as fred
import Fred
import numpy as np

m = 50
sep = 10.

p1 = np.array((m-1) * [[0, 0]] + [[sep, 0]])
p2 = np.array([[0, -0.1], [-2*sep, -0.1]] + (m-2) * [[0, -0.1]])

curves = fred.Curves()


curves.add(fred.Curve(p1, "Input Curve {}".format(1)))
curves.add(fred.Curve(p2, "Input Curve {}".format(2)))


clustering = fred.two_two_dtw_one_two_median(curves)
clustering_e = fred.two_two_dtw_one_two_median_exact(curves)

Fred.plot_curve(curves, clustering, clustering_e)
