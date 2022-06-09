from .backend import *

config = Config()

def plot_curve(*curves, vertex_markings=True, savename=None, saveextension=None, return_fig=False, legend=True):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    max_compl = 1
    max_dim = 1
    fig = plt.figure()
    ax = None
    for curve in curves:
        if isinstance(curve, backend.Curve):
            max_compl = max(max_compl, curve.complexity)
            max_dim = max(max_dim, curve.dimensions)
        elif isinstance(curve, backend.Curves):
            for curv in curve:
                max_compl = max(max_compl, curv.complexity)
                max_dim = max(max_dim, curv.dimensions)
        elif isinstance(curve, backend.Clustering_Result):
            for curv in curve:
                max_compl = max(max_compl, curv.complexity)
                max_dim = max(max_dim, curv.dimensions)
    if max_dim >= 3:
        ax = fig.gca(projection='3d')
    else:
        ax = fig.gca()
    for curve in curves:
        if isinstance(curve, backend.Curve):
            if curve.dimensions >= 3:
                p = ax.plot(curve.values[:, 0], curve.values[:, 1], curve.values[:, 2], linestyle='--', marker='o', label = curve.name, markersize = 7, markevery = curve.complexity)
                if vertex_markings:
                    ax.plot(curve.values[1:, 0], curve.values[1:, 1], curve.values[:, 2], linestyle='', marker='x', label = None, color = p[0].get_color(), markersize = 7)
            elif curve.dimensions == 2:
                p = ax.plot(curve.values[:, 0], curve.values[:, 1], linestyle='--', marker='o', label = curve.name, markersize = 7, markevery = curve.complexity)
                if vertex_markings:
                    ax.plot(curve.values[1:, 0], curve.values[1:, 1], linestyle='', marker='x', label = None, color = p[0].get_color(), markersize = 7)
            else:
                p = ax.plot([i * max_compl / len(curve) for i in range(len(curve))], curve.values, linestyle='--', marker='o', label = curve.name, markersize = 7, markevery = curve.complexity)
                if vertex_markings:
                    ax.plot([i * max_compl / len(curve) for i in range(1, len(curve))], curve.values[1:], linestyle='', marker='x', label = None, color = p[0].get_color(), markersize = 7)
        elif isinstance(curve, backend.Curves):
            for curv in curve:
                if curv.dimensions >= 3:
                    p = ax.plot(curv.values[:, 0], curv.values[:, 1], curv.values[:, 2], linestyle='--', marker='o', label = curv.name, markersize = 7, markevery = curv.complexity)
                    if vertex_markings:
                        ax.plot(curv.values[1:, 0], curv.values[1:, 1], curv.values[:, 2], linestyle='', marker='x', label = None, color = p[0].get_color(), markersize = 7)
                elif curv.dimensions == 2:
                    p = plt.plot(curv.values[:, 0], curv.values[:, 1], linestyle='--', marker='o', label = curv.name, markersize = 7, markevery = curv.complexity)
                    if vertex_markings:
                        plt.plot(curv.values[1:, 0], curv.values[1:, 1], linestyle='', marker='x', label = None, color = p[0].get_color(), markersize = 7)
                else:
                    p = plt.plot([i * max_compl / len(curv) for i in range(len(curv))], curv.values, linestyle='--', marker='o', label = curv.name, markersize = 7, markevery = curv.complexity)
                    if vertex_markings:
                        plt.plot([i * max_compl / len(curv) for i in range(1, len(curv))], curv.values[1:], linestyle='', marker='x', label = None, color = p[0].get_color(), markersize = 7)
        elif isinstance(curve, backend.Clustering_Result):
            for curv in curve:
                if curv.dimensions >= 3:
                    p = ax.plot(curv.values[:, 0], curv.values[:, 1], curv.values[:, 2], linestyle='-', marker='o', label = curv.name, markersize = 7, markevery = curv.complexity)
                    if vertex_markings:
                        ax.plot(curv.values[1:, 0], curv.values[1:, 1], curv.values[:, 2], linestyle='', marker='x', label = None, color = p[0].get_color(), markersize = 7)
                elif curv.dimensions == 2:
                    p = plt.plot(curv.values[:, 0], curv.values[:, 1], linestyle='-', marker='o', label = curv.name, markersize = 7, markevery = curv.complexity)
                    if vertex_markings:
                        plt.plot(curv.values[1:, 0], curv.values[1:, 1], linestyle='', marker='x', label = None, color = p[0].get_color(), markersize = 7)
                else:
                    p = plt.plot([i * max_compl / len(curv) for i in range(len(curv))], curv.values, linestyle='-', marker='o', label = curv.name, markersize = 7, markevery = curv.complexity)
                    if vertex_markings:
                        plt.plot([i * max_compl / len(curv) for i in range(1, len(curv))], curv.values[1:], linestyle='', marker='x', label = None, color = p[0].get_color(), markersize = 7)
    if legend:
        ax.legend(title='Curve names:')
    ax.set_title('Fred Curves')
    if savename is not None:
        plt.savefig("{}.{}".format(savename, saveextension), dpi=150)
        plt.close()
    elif return_fig:
        return fig
    else:
        plt.show()
        plt.close()
    
def plot_clustering(clustering_result, curves, vertex_markings=True, savename=None, saveextension=None, return_fig=False, legend=True):
    if not (isinstance(clustering_result, backend.Clustering_Result) and isinstance(curves, backend.Curves)):
        print("Check parameters!")
        return
    if len(clustering_result.assignment) < 1:
        print("compute_assignment was not called! calling now")
        clustering_result.compute_assignment(curves)
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    colors = list(mcolors.BASE_COLORS)
    if len(clustering_result) > len(colors):
        colors = list(mcolors.TABLEAU_COLORS)
    if len(clustering_result) > len(colors):
        colors = list(mcolors.mcolors.CSS4_COLORS)
    max_compl = 1
    max_dim = 1
    fig = plt.figure()
    ax = None
    for curve in curves:
        max_compl = max(max_compl, curve.complexity)
        max_dim = max(max_dim, curve.dimensions)
    if max_dim >= 3:
        ax = fig.gca(projection='3d')
    else:
        ax = fig.gca()
    for i, curve in enumerate(clustering_result):
        if curve.dimensions >= 3:
            p = ax.plot(curve.values[:, 0], curve.values[:, 1], curve.values[:, 2], linestyle='-', marker='o', label = curve.name, color = colors[i], markersize = 7, markevery = curve.complexity)
            if vertex_markings:
                ax.plot(curve.values[1:, 0], curve.values[1:, 1], curve.values[:, 2], linestyle='', marker='x', label = None, color = p[0].get_color(), markersize = 7)
        elif curve.dimensions == 2:
            p = ax.plot(curve.values[:, 0], curve.values[:, 1], linestyle='-', marker='o', label = curve.name, color = colors[i], markersize = 7, markevery = curve.complexity)
            if vertex_markings:
                ax.plot(curve.values[1:, 0], curve.values[1:, 1], linestyle='', marker='x', label = None, color = p[0].get_color(), markersize = 7)
        else:
            p = ax.plot([i * max_compl / len(curve) for i in range(len(curve))], curve.values, linestyle='-', marker='o', label = curve.name, color = colors[i], markersize = 7, markevery = curve.complexity)
            if vertex_markings:
                ax.plot([i * max_compl / len(curve) for i in range(1, len(curve))], curve.values[1:], linestyle='', marker='x', label = None, color = p[0].get_color(), markersize = 7)
    for i in range(len(clustering_result.assignment)):
        for j in range(clustering_result.assignment.count(i)):
            curve = curves[clustering_result.assignment.get(i,j)]
            if curve.dimensions >= 3:
                p = ax.plot(curve.values[:, 0], curve.values[:, 1], curve.values[:, 2], linestyle=':', marker='o', label = curve.name, color = colors[i], markersize = 7, markevery = curve.complexity)
                if vertex_markings:
                    ax.plot(curve.values[1:, 0], curve.values[1:, 1], curve.values[:, 2], linestyle='', marker='x', label = None, color = p[0].get_color(), markersize = 7)
            elif curve.dimensions == 2:
                p = ax.plot(curve.values[:, 0], curve.values[:, 1], linestyle=':', marker='o', label = curve.name, color = colors[i], markersize = 7, markevery = curve.complexity)
                if vertex_markings:
                    ax.plot(curve.values[1:, 0], curve.values[1:, 1], linestyle='', marker='x', label = None, color = p[0].get_color(), markersize = 7)
            else:
                p = ax.plot([i * max_compl / len(curve) for i in range(len(curve))], curve.values, linestyle=':', marker='o', label = curve.name, color = colors[i], markersize = 7, markevery = curve.complexity)
                if vertex_markings:
                    ax.plot([i * max_compl / len(curve) for i in range(1, len(curve))], curve.values[1:], linestyle='', marker='x', label = None, color = p[0].get_color(), markersize = 7)
    if legend:
        ax.legend(title='Curve names:')
    ax.set_title('Fred Clustering')
    if savename is not None:
        plt.savefig("{}.{}".format(savename, saveextension), dpi=150)
        plt.close()
    elif return_fig:
        return fig
    else:
        plt.show()
        plt.close()
