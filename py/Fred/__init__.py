from . import backend

import numpy as np
import matplotlib.pyplot as plt


def plot_curve(*curves):
    for curve in curves:
        points = list()
        for i in range(0, len(curve)):
            points.append(curve[i].values)
        points = np.array(points)
        plt.plot(points[:, 0], points[:, 1])
    plt.show()
