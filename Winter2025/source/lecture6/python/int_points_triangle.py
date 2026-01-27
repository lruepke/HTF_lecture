"""
DEPRECATED: This module is deprecated. Use fem_integration.py instead.

This module has been replaced by fem_integration.py which provides:
- Support for many more quadrature rules (1, 3, 4, 6, 7, 12, 25, 77 points)
- Comprehensive documentation
- Accuracy verification

Migration:
    OLD: from int_points_triangle import int_points_triangle
    NEW: from fem_integration import ip_triangle
         OR from fem_integration import int_points_triangle (legacy wrapper)
"""

import numpy as np
import warnings

def int_points_triangle(nip) :
    """DEPRECATED: Use fem_integration.ip_triangle instead."""
    warnings.warn(
        "int_points_triangle.py is deprecated. Use fem_integration.py instead:\n"
        "  from fem_integration import ip_triangle",
        DeprecationWarning,
        stacklevel=2
    )

    if nip==3 :
        gauss   = np.array([[ 1/6, 2/3, 1/6], [1/6, 1/6, 2/3]]).T.copy()
        weights = np.array([1/6, 1/6, 1/6])
    elif nip==6 :
        gauss = np.zeros((nip,2))
        g1 = (8-np.sqrt(10) + np.sqrt(38-44*np.sqrt(2/5)))/18
        g2 = (8-np.sqrt(10) - np.sqrt(38-44*np.sqrt(2/5)))/18
        gauss[0,0] = 1-2*g1                 #0.108103018168070;
        gauss[0,1] = g1                      #0.445948490915965;
        gauss[1,0] = g1                      #0.445948490915965;
        gauss[1,1] = 1-2*g1                  #0.108103018168070;
        gauss[2,0] = g1                      #0.445948490915965;
        gauss[2,1] = g1                      #0.445948490915965;
        gauss[3,0] = 1-2*g2                  #0.816847572980459;
        gauss[3,1] = g2                      #0.091576213509771;
        gauss[4,0] = g2                      #0.091576213509771;
        gauss[4,1] = 1-2*g2                  #0.816847572980459;
        gauss[5,0] = g2                        #0.091576213509771;
        gauss[5,1] = g2                        #0.091576213509771;

        weights = np.zeros(nip)
        w1 = (620+np.sqrt(213125-53320*np.sqrt(10)))/3720
        w2 = (620-np.sqrt(213125-53320*np.sqrt(10)))/3720
        weights[0] =  w1;                       #0.223381589678011;
        weights[1] =  w1;                       #0.223381589678011;
        weights[2] =  w1;                       #0.223381589678011;
        weights[3] =  w2;                       #0.109951743655322;
        weights[4] =  w2;                       #0.109951743655322;
        weights[5] =  w2;                       #0.109951743655322;
        weights     = 0.5*weights
    else :
        print("not implemented")

    return gauss,weights