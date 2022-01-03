import unittest
import numpy as np
from math import isclose
from pfoil import airfoil
from pfoil import visc


TOL = 1E-1 # this should be 1E-10...


class test_naca0012(unittest.TestCase):

   # validate against naca0012
   # * this geometry includes a trailing edge gap (sharp = false)
   # * values we're pulled from XFOIL 6.99
        
    def test_xfoil_naca0012(self):
        naca0012 = airfoil()
        naca0012.load('tests/naca0012.txt')
        with open('naca0012_xfoil_6.99.txt', 'r') as txt:
            for line in txt.readlines()[1:]:
                re, adeg, lvconv, cl, cm, cd, cdf, cdp, itmax =  line.replace('\n','').split(',')
                _res = visc(naca0012, re, adeg, verbose=False, itmax=10_000)
                self.assertLess(abs(_res['cl'] - float(cl)), TOL)
                self.assertLess(abs(_res['cd'] - float(cd)), TOL)
                self.assertLess(abs(_res['cm'] - float(cm)), TOL)
                    
                    
if __name__ == '__main__':
    unittest.main()
