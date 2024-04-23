import numpy as np

def shapes_tri(xi, eta):
    
    #shape functions
    eta2 = xi
    eta3 = eta
    eta1 = 1-eta2-eta3

    N1 = eta1
    N2 = eta2
    N3 = eta3

    N = np.array([N1, N2, N3])
    
    # and their derivatives
    dNds = np.zeros((2,3))
 
    dNds[0,0]   =  -1 #derivative 
    dNds[1,0]   =  -1 #derivative 

    #derivatives of second shape function with local coordinates
    dNds[0,1]   =  1
    dNds[1,1]   =  0

    #derivatives of third shape function with local coordinates
    dNds[0,2]   =  0
    dNds[1,2]   =  1
    
    return N, dNds