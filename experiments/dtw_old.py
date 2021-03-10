import Fred.backend as fred
import Fred
import numpy as np


c1 = fred.Curve(np.array([[0.0, -1.0], [-1.0, 3.0], [1.0, 0.1], [6.1, 0.1], [6.2, 0.1], [6.3, 0.1], [6.4, 0.1], [6.5, 0.1], [7.0, 0.1], [8.0, 0.1], [8.3, 0.1], [8.4, 0.1], [8.5, 0.1], [8.6, 0.1], [8.7, 0.1], [8.8, 0.1], [8.9, 0.1], [10, -1.0]]), "Eingabe 1")
c2 = fred.Curve(np.array([[0.0, -1.0], [1.0, -0.1], [1.1, -0.1], [1.2, -0.1], [1.3, -0.1], [1.4, -0.1], [1.5, -0.1], [1.6, -0.1], [2.0, -0.1], [3.0, -0.1], [4.0, -0.1], [5.0, -0.1], [6.0, -0.1], [7.0, -0.1], [8.0, -0.1], [11.0, 3.0], [10.0, -1.0]]), "Eingabe 2")

curves = fred.Curves()
curves.add(c1)
curves.add(c2)


clustering = fred.two_two_dtw_one_two_median(curves)
clustering_e = fred.two_two_dtw_one_two_median_exact(curves)

print(clustering.value/clustering_e.value)

Fred.plot_curve(curves, clustering, clustering_e)
