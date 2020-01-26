import numpy as np
import unittest
import Fred

eps = 0.00001

class TestContinuousFrechet(unittest.TestCase):

    def test_zigzag1d(self):
        a = Fred.Curve(np.array([[0.0], [1.0], [0.0], [1.0]]))
        b = Fred.Curve(np.array([[0.0],[0.75], [0.25], [1.0]]))
        c = Fred.Curve(np.array([[0.0],[1.0]]))
        self.assertEqual(Fred.continuous_frechet(a, b, eps).value, 0.25)
        self.assertEqual(Fred.continuous_frechet(a, c, eps).value, 0.5)
        
    def test_longsegment(self):
        a = Fred.Curve(np.array([[0.0, 0.0],[500.0e3, 0.0], [1.0e6, 0.0]]))
        b = Fred.Curve(np.array([[0.0, 0.0],[1.0e6, 0.0]]))
        self.assertEqual(Fred.continuous_frechet(a, b, eps).value, 0.0)

class TestDiscreteFrechet(unittest.TestCase):
    
    def test_zigzag1d(self):
        a = Fred.Curve(np.array([[0.0], [1.0], [0.0], [1.0]]))
        b = Fred.Curve(np.array([[0.0],[0.75], [0.25], [1.0]]))
        c = Fred.Curve(np.array([[0.0],[1.0]]))
        self.assertEqual(Fred.discrete_frechet(a, b).value, 0.25)
        self.assertEqual(Fred.discrete_frechet(a, c).value, 1.0)
        
    def test_longsegment(self):
        a = Fred.Curve(np.array([[0.0, 0.0],[500.0e3, 0.0], [1.0e6, 0.0]]))
        b = Fred.Curve(np.array([[0.0, 0.0],[1.0e6, 0.0]]))
        self.assertEqual(Fred.discrete_frechet(a, b).value, 500000.0)

if __name__ == '__main__':
    unittest.main()
