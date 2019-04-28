#!/usr/bin/python
from __future__ import print_function
import numpy as np

def test_save():
    rext = 1.1
    Ns = 4
    x, y = np.mgrid[ -rext:rext:1j*Ns, -rext:rext:1j*Ns ]
    
    print( "shape ", x.shape, y.shape )
    
    np.save("tmp", (x,y) ) 
    
    newx,newy = np.load("tmp.npy")
    
    print( "recover shape ", newx.shape, newy.shape)
    assert( newx.shape == x.shape )
    assert( newy.shape == y.shape )
    
    print( "orig x\n",x )
    print( "recovered x\n",newx )
    assert( newx[3,1] == x[3,1] )
    assert( newx[2,0] == x[2,0] )
    
    print( "orig y\n",y )
    print( "recovered x\n",newy )
    assert( newy[3,1] == y[3,1] )
    assert( newy[2,3] == y[2,3] )
    
    return

# ==================================
if __name__ == '__main__':
    test_save()
