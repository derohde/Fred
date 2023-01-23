"""
Copyright 2022 - 2023 Alexander Neuhaus, Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
import numpy as np
import numba
import math
from scipy import optimize
import random

def rejection_sampling(dim, radius, center, n):
    points = []
    if radius == 0:
        points.append(center)
    else:
        for i in range(n):
            while True:
                p = np.random.uniform(low=0.0, high=1.0, size=dim) * 2 - 1
            
                if np.linalg.norm(p) <= 1:
                    break
            p *= radius
            p += center
            points.append(p)
    return np.array(points)

def unit_vector(vector):
    if np.linalg.norm(vector) == 0:
        return vector
    return vector / np.linalg.norm(vector)

def angle(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)

    return np.pi - abs(np.arccos(np.dot(v1_u, v2_u)))

def f(x, *args):
    if len(x) < 3:
        return -np.pi
    samples = args[0]
    start = args[1]
    end =  args[2]
    y = [int(xs) for xs in x]
    curve = list()
    vectors = list()
    angles = list()
    for i in range(start, end):
        curve.append(samples[i][y[i-start]])
    for i in range(len(curve)-1):
        vectors.append(curve[i+1] - curve[i])
    for i in range(len(vectors)-1):
        angles.append(-angle(vectors[i],vectors[i+1]))
    return max(angles)

def compute_stabber(samples, start, end):
    lb = list()
    ub = list()
    for i in range(start, end):
        lb.append(0)
        ub.append(len(samples[i])-0.001)

    res = optimize.dual_annealing(f, args=(samples, start, end), bounds=list(zip(lb,ub)))
    xf = res.x
    curve = list()
    for i in range(start, end):
        curve.append(samples[i][int(xf[i-start])])
    return(np.array(curve))

def outer_tangents(r1, c1, r2, c2):
    angle1 = math.atan2(c2[1] - c1[1], c2[0] - c1[0])
    if(c1[0] == c2[0] and c1[1] == c2[1]):
        angle2 = math.pi/2
    else:
        angle2 = math.acos(np.clip((r1 - r2) / np.linalg.norm(c1 - c2), -1, 1))

    t1StartX = c1[0] + math.cos(angle1 + angle2) * r1
    t1StartY = c1[1] + math.sin(angle1 + angle2) * r1
    t1EndX = c2[0] + math.cos(angle1 + angle2) * r2
    t1EndY = c2[1] + math.sin(angle1 + angle2) * r2

    t2StartX = c1[0] + np.cos(angle1 - angle2) * r1
    t2StartY = c1[1] + np.sin(angle1 - angle2) * r1
    t2EndX = c2[0] + np.cos(angle1 - angle2) * r2
    t2EndY = c2[1] + np.sin(angle1 - angle2) * r2

    t1 = [[t1StartX, t1StartY], [t1EndX, t1EndY]]
    t2 = [[t2StartX, t2StartY], [t2EndX, t2EndY]]
    return t1, t2

def convert_coordinates(p0, p1, p2):
    if len(p0) < 3:
        l = len(p0)
        p0t = np.zeros(3)
        p1t = np.zeros(3)
        p2t = np.zeros(3)
        for i in range(l):
            p0t[i] = p0[i]
            p1t[i] = p1[i]
            p2t[i] = p2[i]
        p0, p1, p2 = p0t, p1t, p2t
        
    x = unit_vector(p1 - p0)
    z = np.cross(x,(p2 - p0))
    y = np.cross(z, x)
    o = p0
    
    #projection of p0 to the origin
    x0 = 0
    y0 = 0
    
    #projection of p1 to d(p1,p0),0
    x1 = np.linalg.norm(p1 - o)
    y1 = 0
    
    #projection of p2
    x2 = np.dot((p2 - o), x)
    y2 = np.dot((p2 - o), y)
    
    return np.array((x0,y0)), np.array((round(x1, 12), round(y1, 12))), np.array((round(x2, 12), round(y2, 12)))

def check_containment(c1, r1, c2, r2, p):
    p1, p2, p3 = convert_coordinates(c1, c2, p)
    t1, t2 = outer_tangents(r1, p1, r2, p2)
    t1A = t1[0]
    t1B = t1[1]
    t2A = t2[0]
    t2B = t2[1]
    
    p3p = np.array((p3[0], 0))
    
    if((np.linalg.norm(p - c1) <= r1) or (np.linalg.norm(p - c2) <= r2)):
        return True
    else:
        sign1 = (p3[0] - t1A[0]) * (t1B[1] - t1A[1]) - (p3[1] - t1A[1]) * (t1B[0] - t1A[0])
        sign2 = (p3[0] - t2A[0]) * (t2B[1] - t2A[1]) - (p3[1] - t2A[1]) * (t2B[0] - t2A[0])
        if((sign1 >= 0 and sign2 <= 0)
           and (np.linalg.norm(p3p - p1) <= np.linalg.norm(p1 - p2))
           and (np.linalg.norm(p3p - p2) <= np.linalg.norm(p1 - p2))):
            return True
        return False

def is_stabbable(balls, old_samples, new_samples, start, end):
    for j in range(start + 1, end):
        for i in range(start, j):
            c1 = balls[i][0]
            r1 = balls[i][1]
            c2 = balls[end][0]
            r2 = balls[end][1]
            points = [p for p in old_samples[j] if check_containment(c1.copy(), r1, c2.copy(), r2, p.copy())]
            new_samples[j] = np.array(points)
            if len(points) == 0:
                return False
    return True

def stabbing_path(balls, epsilon=0.1, n_samples=None):
    dim = len(balls[0][0])
    if n_samples is None:
        n_samples = int(10 * 1/epsilon * np.log(len(balls)))
    old_samples = list()
    new_samples = list()
    segments = list()
    for i in range(len(balls)):
        new_samples.append(rejection_sampling(dim, balls[i][1], balls[i][0], n_samples))
        old_samples.append(new_samples[i])
    start = 0;
    end = 0;
    stabbable = True
    while end < len(balls):
        if stabbable:            
            old_samples = new_samples.copy()
            stabbable = is_stabbable(balls, old_samples, new_samples, start, end)
            end += 1
        else:
            end -= 1
            tmp = compute_stabber(old_samples, start, end-1)
            segments.append((tmp[0],tmp[len(tmp)-1]))
            new_samples = old_samples.copy()
            start = end
            stabbable = True
    if not stabbable:
        tmp = compute_stabber(old_samples, start, end - 1)
        segments.append((tmp[0],tmp[len(tmp) - 1]))
        segments.append(np.array([old_samples[-1][0]]))
    else:
        tmp = compute_stabber(old_samples, start, end)
        segments.append((tmp[0],tmp[len(tmp)-1]))
        old_samples = new_samples
    return np.concatenate(segments), old_samples
