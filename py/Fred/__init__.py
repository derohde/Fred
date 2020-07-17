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
                plt.plot(points[:, 0], points[:, 1], label = curve.name)
            else:
                plt.plot(points, label = curve.name)
        elif isinstance(curve, backend.Curves):
            for i in range(0, len(curve)):
                points = list()
                for j in range(0, len(curve[i])):
                    points.append(curve[i][j].values)
                points = np.array(points)
                if curve[i].dimensions >= 2:
                    plt.plot(points[:, 0], points[:, 1], label = curve[i].name)
                else:
                    plt.plot(points, label = curve[i].name)
    plt.legend(title='Curve names:')
    plt.title('Fred Curves')
    plt.show()
