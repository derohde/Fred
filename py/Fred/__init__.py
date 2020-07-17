from . import backend

import numpy as np
import matplotlib.pyplot as plt


def plot_curve(*curves):
    for curve in curves:
        if isinstance(curve, backend.Curve):
            points = list()
            for i in range(0, len(curve)):
                points.append(curve[i].values)
            points = np.array(points)
            if curve.dimensions >= 2:
                plt.plot(points[:, 0], points[:, 1])
            else:
                plt.plot(points)
        elif isinstance(curve, backend.Curves):
            for i in range(0, len(curve)):
                points = list()
                for j in range(0, len(curve[i])):
                    points.append(curve[i][j].values)
                points = np.array(points)
                if curve.dimensions >= 2:
                    plt.plot(points[:, 0], points[:, 1])
                else:
                    plt.plot(points)
    plt.show()
