from .backend import *

config = Config()

def plot_curve(*curves, savename=None, saveextension=None):
    import matplotlib.pyplot as plt
    max_compl = 1
    for curve in curves:
        if isinstance(curve, backend.Curve):
            max_compl = max(max_compl, curve.complexity)
        elif isinstance(curve, backend.Curves):
            for curv in curve:
                max_compl = max(max_compl, curv.complexity)
        elif isinstance(curve, backend.Clustering_Result):
            for curv in curve:
                max_compl = max(max_compl, curv.complexity)
    for curve in curves:
        if isinstance(curve, backend.Curve):
            if curve.dimensions >= 2:
                p = plt.plot(curve.values[:, 0], curve.values[:, 1], '--o', label = curve.name, markersize = 7, markevery = curve.complexity)
                plt.plot(curve.values[1:, 0], curve.values[1:, 1], 'x', label = None, color = p[0].get_color(), markersize = 7)
            else:
                p = plt.plot([i * max_compl / len(curve) for i in range(len(curve))], curve.values, '--o', label = curve.name, markersize = 7, markevery = curve.complexity)
                plt.plot([i * max_compl / len(curve) for i in range(1, len(curve))], curve.values[1:], 'x', label = None, color = p[0].get_color(), markersize = 7)
        elif isinstance(curve, backend.Curves):
            for curv in curve:
                if curv.dimensions >= 2:
                    p = plt.plot(curv.values[:, 0], curv.values[:, 1], '--o', label = curv.name, markersize = 7, markevery = curv.complexity)
                    plt.plot(curv.values[1:, 0], curv.values[1:, 1], 'x', label = None, color = p[0].get_color(), markersize = 7)
                else:
                    p = plt.plot([i * max_compl / len(curv) for i in range(len(curv))], curv.values, '--o', label = curv.name, markersize = 7, markevery = curv.complexity)
                    plt.plot([i * max_compl / len(curv) for i in range(1, len(curv))], curv.values[1:], 'x', label = None, color = p[0].get_color(), markersize = 7)
        elif isinstance(curve, backend.Clustering_Result):
            for curv in curve:
                if curv.dimensions >= 2:
                    p = plt.plot(curv.values[:, 0], curv.values[:, 1], '-o', label = curv.name, markersize = 7, markevery = curv.complexity)
                    plt.plot(curv.values[1:, 0], curv.values[1:, 1], 'x', label = None, color = p[0].get_color(), markersize = 7)
                else:
                    p = plt.plot([i * max_compl / len(curv) for i in range(len(curv))], curv.values, '-o', label = curv.name, markersize = 7, markevery = curv.complexity)
                    plt.plot([i * max_compl / len(curv) for i in range(1, len(curv))], curv.values[1:], 'x', label = None, color = p[0].get_color(), markersize = 7)
    plt.legend(title='Curve names:')
    plt.title('Fred Curves')
    if savename is None:
        plt.show()
    else:
        plt.savefig("{}.{}".format(savename, saveextension), dpi=150)
    plt.close()
