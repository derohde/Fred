import Fred.backend as fred
import Fred
import numpy as np


seperation = 100
ninb = 20
clusterradius = 5
clusternumber = 20
outlierradius = 200
clusteroutliers = 2
ncurves = 3

curves = fred.Curves()

for i in range(ncurves):
    
    c1 = np.random.normal((0, 0), clusterradius, (clusternumber, 2))
    c2 = np.random.normal((seperation, 0), clusterradius, (clusternumber, 2))
    
    inb = None
    
    for j in range(1, ninb + 1):
        if j == 1:
            inb = np.random.normal(0 + j * seperation / (ninb+2), 1, (1,2))
        else:
            inb = np.concatenate([inb, np.random.normal((0 + j * seperation / (ninb+2), 1), 1, (1,2))])
    
    k = None
    
    for j in range(clusteroutliers):
        kk = int(np.random.rand() * clusternumber)
        while kk == k:
            kk = int(np.random.rand() * clusternumber)
        k = kk
        randomdirection = np.random.normal(0, 1, (1, 2)).flatten()
        randomdirection /= np.linalg.norm(randomdirection)
        c1[k] += outlierradius * randomdirection
    
    k = None
    
    for j in range(clusteroutliers):
        kk = int(np.random.rand() * clusternumber)
        while kk == k:
            kk = int(np.random.rand() * clusternumber)
        k = kk
        randomdirection = np.random.normal(0, 1, (1, 2)).flatten()
        randomdirection /= np.linalg.norm(randomdirection)
        c2[k] += outlierradius * randomdirection
        
    curves.add(fred.Curve(np.concatenate([c1, inb, c2]), "Input Curve {}".format(i+1)))

clustering = fred.two_two_dtw_one_two_median(curves)
clustering_e = fred.two_two_dtw_one_two_median_exact(curves)

Fred.plot_curve(curves, clustering, clustering_e)
    
