import numpy as np
import Fred

eps = 0.00001

a = Fred.Curve(np.array([[0.0], [1.0], [0.0], [1.0]]))
b = Fred.Curve(np.array([[0.0],[0.75], [0.25], [1.0]]))
c = Fred.Curve(np.array([[0.0],[1.0]]))

print(Fred.continuous_frechet(a, b, eps).value) # 0.125
print(Fred.continuous_frechet(a, c, eps).value) # 0.5
print(Fred.discrete_frechet(a, b).value) # 0.125
print(Fred.discrete_frechet(a, c).value) # 0.5

