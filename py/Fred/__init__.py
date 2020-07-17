from . import backend

import numpy as np
import matplotlib.pyplot as plt


def plot_curve(*curves, save=None):
    for curve in curves:
        if isinstance(curve, backend.Curve):
            if curve.dimensions >= 2:
                plt.plot(curve.values[:, 0], curve.values[:, 1], label = curve.name)
            else:
                plt.plot(curve.values, label = curve.name)
        elif isinstance(curve, backend.Curves):
            for curv in curve:
                if curv.dimensions >= 2:
                    plt.plot(curv.values[:, 0], curv.values[:, 1], label = curv.name)
                else:
                    plt.plot(curv.values, label = curv.name)
        elif isinstance(curve, backend.Clustering_Result):
            for center in curve:
                if center.dimensions >= 2:
                    plt.plot(center.values[:, 0], center.values[:, 1], label = center.name)
                else:
                    plt.plot(center.values, label = center.name)
                
    plt.legend(title='Curve names:')
    plt.title('Fred Curves')
    if save is None:
        plt.show()
    else:
        plt.savefig("{}.svg".format(save), dpi=150)
