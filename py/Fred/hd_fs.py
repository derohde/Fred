from . import backend

from cvxopt import matrix, solvers
solvers.options['show_progress'] = False
import numpy as np

def curve_within_radii(T, radii):
    l = len(T)
    start_point = _intersection_point(T, radii, [1] * l, [0] * l)
    if not start_point[0]:
        return (False, )
    end_point = _intersection_point(T, radii, [T[i].points - 1 for i in range(0, l)], [1] * l)
    if not end_point[0]:
        return (False, )
    borders = [([1] * len(T), start_point[1]), ([T[i].points for i in range(0, l)], end_point[1])]
    _compute_borders(T, radii, np.array([1] * l, dtype=np.uint64), borders)
    borders = sorted(borders, key=lambda border: border[0])
    solution = list()
    solution.append(borders[-1])
    alpha = [T[i].points for i in range(0, l)]
    for i in range(len(borders)-2, -1, -1):
        diff = np.array(alpha) - np.array(borders[i][0])
        bools = any(diff > 1) or any(np.array(borders[i][0]) > np.array(alpha)) or all(np.array(borders[i][0]) == np.array(alpha))
        if not bools:
            solution.append(borders[i])
            alpha = borders[i][0]
    if solution[-1] != borders[0]:
        return (False, )
    solution = reversed(solution)
    return (True, backend.Curve(np.array([sol[1] for sol in solution], dtype="double")))
    
def _compute_borders(T, radii, alpha, borders):
    for j in range(0, len(T)):
        if alpha[j] + 1 < T[j].points:
            border = _cellborder(T, radii, alpha, j, 1)
            if border[0]:
                mask = np.zeros((len(T), ))
                mask[j] = 1
                borders.append(([alpha[i] + border[1][0][i] for i in range(0, len(alpha))], border[1][1]))
                _compute_borders(T, radii, alpha + mask, borders)

def _intersection_point(T, radii, alpha, values):
    lines = list()
    for i in range(0, len(T)):
        lines.append([(1-values[i])*T[i][alpha[i]-1][k] + values[i]*T[i][alpha[i]][k] for k in range(0, T[i].dimensions)])
    
    c = matrix([0.] * T[0].dimensions)
    G = []
    for i in range(0, T[0].dimensions):
        m = len(T) * ([0.] + i * [0.] + [-1] + (T[0].dimensions - i - 1) * [0.])
        G.append(m)
        
    G = matrix(G)
    h = []
    for i in range(0, len(T)):
        h += [radii[i]] + [-lines[i][j] for j in range(0, T[i].dimensions)]
    h = matrix(h)
    dims = {'l': 0, 'q': [1+T[0].dimensions] * len(T), 's': [0]}
    sol = solvers.conelp(c, G, h, dims)
    return sol['status'] == 'optimal', np.array([s for s in sol['x']]) if sol['x'] is not None else None

def _cellborder(T, radii, alpha, j, value):
    lines = list()
    for i in range(0, len(T)):
        if (i == j): 
            lines.append(([(1-value)*T[i][int(alpha[i]-1)][k] + value*T[i][int(alpha[i])][k] for k in range(0, T[i].dimensions)],))
        else:
            lines.append(([T[i][int(alpha[i]-1)][k] - T[i][int(alpha[i])][k] for k in range(0, T[i].dimensions)], 
                          [T[i][int(alpha[i]-1)][k] for k in range(0, T[i].dimensions)]))
        
    dims = {'l': 2 * (len(T)-1), 'q': [1+T[0].dimensions] * len(T), 's': [0]}

    c = matrix([0.] * (len(T) - 1) + [0.] * T[0].dimensions)
    G = list()
    
    for i in range(0, len(T)):
        if i == j:
            pass
        elif i < j:
            m = i * [0.] + [1.] + (len(T) - i - 2) * [0.]
            m += i * [0.] + [-1.] + (len(T) - i - 2) * [0.]
            m += (i * (T[i].dimensions + 1)) * [0.] + [0.] + [-lines[i][0][k] for k in range(0, T[i].dimensions)] + ((len(T) - i - 1) * (T[i].dimensions + 1)) * [0.]
            G.append(m)
        else:
            m = (i - 1) * [0.] + [1.] + (len(T) - i - 1) * [0.]
            m += (i - 1) * [0.] + [-1.] + (len(T) - i - 1) * [0.]
            m += (i * (T[i].dimensions + 1)) * [0.] + [0.] + [-lines[i][0][k] for k in range(0, T[i].dimensions)] + ((len(T) - i - 1) * (T[i].dimensions + 1)) * [0.]
            G.append(m)
        
    for i in range(0, T[0].dimensions):
        m = (2 * (len(T) - 1)) * [0.]
        m += len(T) * ([0.] + i * [0.] + [-1] + (T[0].dimensions - i - 1) * [0.])
        G.append(m)
        
    G = matrix(G)
    
    h = [1.] * (len(T) - 1)
    h += [0.] * (len(T) - 1)
    for i in range(0, len(T)):
        h += [radii[i]]
        if i == j:
            h += [-lines[i][0][k] for k in range(0, T[i].dimensions)]
        else:
            h += [-lines[i][1][k] for k in range(0, T[i].dimensions)]
    h = matrix(h)
    try:
        sol = solvers.conelp(c, G, h, dims)
    except Exception as e:
        print(e)
    if sol['x'] is not None:
        x = sol['x']
        y = np.empty((len(T),))
        for i in range(0, len(T)):
            if i < j:
                y[i] = x[i]
            elif i == j:
                y[j] = value
            else:
                y[i] = x[i-1]
    return sol['status'] == 'optimal', (y, np.array([s for s in sol['x'][len(T)-1:]])) if sol['x'] is not None else None
