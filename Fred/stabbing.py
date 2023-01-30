"""
Copyright 2022 - 2023 Alexander Neuhaus, Dennis Rohde

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
import numpy as np
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
    if((c1 == c2).all()):
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
    
def convert_coordinates(p1, p2, p3):
    
    def angle(u, v):
        return np.arccos(np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)))

    def rotate(p, i, ang):
        q = p.copy()
        q[i] = p[i] * np.cos(ang) - p[i + 1] * np.sin(ang)
        q[i + 1] = p[i] * np.sin(ang) + p[i + 1] * np.cos(ang)
        return q
        
    def ax_ang(p, i):
        u = np.zeros_like(p)
        u[i] = 1
        pp = np.zeros_like(p)
        pp[i] = p[i]
        pp[i + 1] = p[i + 1]
        if pp[i + 1] >= 0:
            return np.pi / 2 - angle(pp, u)
        else:
            return np.pi / 2 + angle(pp, u)
        
    def transform(y, z):
        for i in range(0, len(y) - 1):
            ang = ax_ang(y, i)
            y = rotate(y, i, ang)
            z = rotate(z, i, ang)
        return y, z
        
    def transform2(z):
        for i in range(0, len(z) - 2):
            z = rotate(z, i, ax_ang(z,i))
        return z
        
    x = np.zeros_like(p1)
    y = p2 - p1
    z = p3 - p1
    
    y, z = transform(y, z)
    z = transform2(z)
    
    return x[-2:], y[-2:], z[-2:]

def check_containment(c1, r1, c2, r2, p):
    
    def collinear(p1, p2, p3):
        x = np.linalg.norm(p1 - p2)
        y = np.linalg.norm(p2 - p3)
        z = np.linalg.norm(p3 - p1)
        d = [x, y, z]
        m = np.argmax(d)
        dd = [d[i] for i in range(3) if i != m]
        if round(d[m], 12) == round(sum(dd), 12):
            return True
        else:
            return False
        
    if collinear(c1, c2, p):
        if((np.linalg.norm(p - c1) <= r1) or (np.linalg.norm(p - c2) <= r2)):
            return True
        else:
            if((np.linalg.norm(p - c1) <= np.linalg.norm(c1 - c2))
               and (np.linalg.norm(p - c2) <= np.linalg.norm(c1 - c2))):
                return True
            return False
    
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

def stabbing_path(balls, epsilon=0.5, n_samples=None):
    dim = len(balls[0][0])
    if n_samples is None:
        n_samples = int(100 * 1/epsilon * np.log(len(balls)))
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
