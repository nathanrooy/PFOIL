import numpy as np
from ctypes import c_float
from xfoil import pfoil_visc


def aread(airfoil_path):    

    i, ibx = 0, 1480
    x = np.asfortranarray(np.zeros(ibx), dtype=c_float)
    y = np.asfortranarray(np.zeros(ibx), dtype=c_float)
        
    with open(airfoil_path, 'r') as txt:
        for line in txt.readlines():
            try:
                if len(line.split())==2:
                    _x, _y = list(map(float, line.split()))
                    x[i] = _x
                    y[i] = _y
                    i += 1
            except:
                pass
    
    name = airfoil_path.split('/')[-1].split('.')[0].upper()
    return x, y, i, name
    
    
class airfoil():
    def __init__(self):
        pass
    
    
    def load(self, airfoil_path):
        '''Load airfoil from file
        '''
        # check if file exists first...
        self.airfoil_path = airfoil_path
        
        # load airfoil geometry
        self.x, self.y, self.n, self.name = aread(self.airfoil_path)
        
        # compute initial geometric quantities
        
        
    def naca(self, naca_digits):
    
        # compute initial geometric quantities
    
        pass
        
        
    def bezier(self, cps):
    
        # compute initial geometric quantities
    
        pass
        
        

def visc(foil, re, adeg, itmax=20, verbose=True, return_pfs_only=False):
    
    cl, cm, cd, cdf, cdp = pfoil_visc(foil.x, foil.y, float(adeg), float(re), int(foil.n), int(itmax), bool(verbose))
   
   
    return {'cl':cl, 'cm':cm, 'cd':cd, 'cdf':cdf, 'cdp':cdp}
