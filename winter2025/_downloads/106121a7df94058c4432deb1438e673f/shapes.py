import numpy as np

def shapes(xi, eta):
    
    #shape functions
    
    N1 = 0.25*(1-xi)*(1-eta)
    N2 = 0.25*(1+xi)*(1-eta)
    N3 = 0.25*(1+xi)*(1+eta)
    N4 = 0.25*(1-xi)*(1+eta)

    N = np.array([N1, N2, N3, N4])
    
    # and their derivatives
    dNds = np.zeros((2,4))
 
    dNds[0,0]   =  0.25*(-1  + eta) #derivative with xi
    dNds[1,0]   =  0.25*(xi  -  1) #derivative with eta

    #derivatives of second shape function with local coordinates
    dNds[0,1]   =  0.25*(1   - eta)
    dNds[1,1]   =  0.25*(-xi -  1)

    #derivatives of third shape function with local coordinates
    dNds[0,2]   =  0.25*(eta  +  1)
    dNds[1,2]   =  0.25*(xi  +  1)

    #derivatives of fourth shape function with local coordinates
    dNds[0,3]   =  0.25*(-eta -  1)
    dNds[1,3]   =  0.25*(1   - xi)
    
    return N, dNds