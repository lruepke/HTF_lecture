import numpy as np

def shapes_tri(xi, eta, nnodel):
    #shape functions
    eta2 = xi
    eta3 = eta
    eta1 = 1-eta2-eta3

    if nnodel==3 :
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
        
    elif nnodel==6 :

        N    = np.array([eta1*(2*eta1-1),eta2*(2*eta2-1),eta3*(2*eta3-1),4*eta2*eta3,4*eta1*eta3,4*eta1*eta2])       
        dNds = np.array([[1-4*eta1, -1+4*eta2,0, 4*eta3, -4*eta3,4*eta1-4*eta2],[1-4*eta1,0, -1+4*eta3, 4*eta2,4*eta1-4*eta3,-4*eta2]])

    return N, dNds