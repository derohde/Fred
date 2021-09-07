import numpy as np
import unittest
import Fred.backend as fred

class TestContinuousFrechet(unittest.TestCase):

    def test_zigzag1d(self):
        a = fred.Curve([0.0, 1.0, 0.0, 1.0])
        b = fred.Curve([0.0, 0.75, 0.25, 1.0])
        c = fred.Curve([0.0, 1.0])
        self.assertEqual(round(fred.continuous_frechet(a, b).value, 2), 0.25)
        self.assertEqual(round(fred.continuous_frechet(a, c).value, 1), 0.5)
        
    def test_longsegment(self):
        a = fred.Curve([0.0,500.0e3, 1.0e6])
        b = fred.Curve([0.0, 1.0e6])
        self.assertEqual(round(fred.continuous_frechet(a, b).value, 3), 0.0)

class TestDiscreteFrechet(unittest.TestCase):
    
    def test_zigzag1d(self):
        a = fred.Curve([0.0, 1.0, 0.0, 1.0])
        b = fred.Curve([0.0, 0.75, 0.25, 1.0])
        c = fred.Curve([0.0, 1.0])
        self.assertEqual(fred.discrete_frechet(a, b).value, 0.25)
        self.assertEqual(fred.discrete_frechet(a, c).value, 1.0)
        
    def test_longsegment(self):
        a = fred.Curve([0.0,500.0e3, 1.0e6])
        b = fred.Curve([0.0, 1.0e6])
        self.assertEqual(fred.discrete_frechet(a, b).value, 500000.0)
        
class TestDiscreteDynamicTimeWarping(unittest.TestCase):
    
    def test_zigzag1d(self):
        a = fred.Curve([0.0, 1.0, 0.0, 1.0])
        b = fred.Curve([0.0, 0.75, 0.25, 1.0])
        c = fred.Curve([0.0, 1.0])
        self.assertEqual(fred.discrete_dynamic_time_warping(a, b).value, 0.5)
        self.assertEqual(fred.discrete_dynamic_time_warping(a, c).value, 1.0)
        
    def test_longsegment(self):
        a = fred.Curve([0.0,500.0e3, 1.0e6])
        b = fred.Curve([0.0, 1.0e6])
        self.assertEqual(fred.discrete_dynamic_time_warping(a, b).value, 500000.0)

if __name__ == '__main__':
    unittest.main()
