import numpy as np
import argparse
import os
import matplotlib.pyplot as plt
import triangle as tr
from mechanical2d import *

def make_mesh(x_min, x_max, y_min, y_max, c_inc, no_pts,radius, el_ids = (1,100), el_sizes=(0.1, 0.01)):
    
    def make_box():
        i = len(vertices)

        vertices.extend([[x_min,   y_min],
                        [x_max, y_min],
                        [x_max, y_max],
                        [x_min,   y_max]])

        segments.extend([(i+0, i+1),
                        (i+1, i+2),
                        (i+2, i+3),
                        (i+3, i+0)])
        
        regions.append([x_min+0.01*(x_max-x_min), y_min+0.01*(y_max-y_min), el_ids[0],el_sizes[0]])

        segment_markers.extend([101, 102, 103, 104])

    
    def make_inclusion():
        theta        = np.linspace(0,2*np.pi,no_pts, endpoint=False)
        xx           = np.cos(theta)
        yy           = np.sin(theta)

        i = len(vertices)   
        vertices.extend(np.array([c_inc[0] + radius*xx,c_inc[1] + radius*yy]).T)

        Tmp = np.array([np.arange(i, i+no_pts), np.arange(i+1, i+no_pts+1)]).T
        Tmp[-1,1] = i
        segments.extend(Tmp)
        segment_markers.extend(101*np.ones(len(Tmp)))

        regions.append([c_inc[0], c_inc[1], el_ids[1],el_sizes[1]])

    # arrays to fill in with input
    vertices = []
    segments = []
    regions  = []
    segment_markers = []

    #make geometry
    make_box()
    make_inclusion()

    # create mesh by calling triangle
    A = dict(vertices=vertices, segments=segments, segment_markers=segment_markers, regions=regions)
    B = tr.triangulate(A, 'o2pq33Aa')

    # extract mesh information
    GCOORD = B.get("vertices")
    EL2NOD = B.get("triangles")
    Phases = B.get("triangle_attributes").flatten().astype(int)  #makes a copy
    Node_ids = B.get("vertex_markers").flatten()

    #modify GCOORD and EL2NOD to so that a 7th node is added at the center
    nnod        = GCOORD.shape[0]
    GCOORD      = np.vstack((GCOORD, np.mean(GCOORD[EL2NOD[:,0:3]], axis=1))) #now we made a copy
    new_nodes   = np.arange(nnod, GCOORD.shape[0]).reshape(-1, 1)
    EL2NOD      = np.hstack((EL2NOD, new_nodes))

    #make new Mesh dict
    Mesh = dict(
        GCOORD=GCOORD, 
        EL2NOD=EL2NOD, 
        Phases=Phases,
        nnod=GCOORD.shape[0],
        nel=EL2NOD.shape[0],
        nnodel=EL2NOD.shape[1],
        Node_ids=Node_ids)

    return Mesh

if __name__ == "__main__":
    
    # Define material properties
    D           = np.array([1e-3,  1.])      #Viscosity
    Rho         = np.array([   1.,  2.])      #Density
    G           = np.array([   0.,  0.])      #Gravity

    # Define mesh properties 
    no_pts      =  60
    radius      =  0.2
    c_inc       =  (0,0)
    x_min		= -1
    x_max		=  1
    y_min		= -1
    y_max		=  1
    el_sizes    = (0.01, 0.001)
    el_ids      = (0,1)   

    # and make it
    Mesh = make_mesh(x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, c_inc=c_inc, no_pts=no_pts,radius=radius, el_sizes=el_sizes, el_ids=el_ids)

    #set boundary conditions
    #note that Bc_ind points to degree of freedeom and not node number!
    Bc_ind  = np.where(np.isin(Mesh.get('Node_ids'), [101, 102, 103, 104]))[0]
    Bc_val  = np.hstack((Mesh.get('GCOORD')[Bc_ind,0],  -Mesh.get('GCOORD')[Bc_ind,1]))
    Bc_ind  = np.hstack((2*Bc_ind, 2*Bc_ind+1))

    # and add to Mesh structure
    Mesh['Bc_ind'] = Bc_ind
    Mesh['Bc_val'] = Bc_val

    Materials = dict(D=D, Rho=Rho, G=G)

    Solver = dict(nip=6)

    # now call solver
    Vel = mechanical2d(Mesh=Mesh, Materials=Materials, Solver=Solver  )

    Vel_x = Vel[0::2]
    Vel_y = Vel[1::2]

    #make quiver plot
    plt.figure()
    plt.triplot(Mesh.get('GCOORD')[:,0], Mesh.get('GCOORD')[:,1], Mesh.get('EL2NOD')[:,0:3])
    plt.quiver(Mesh.get('GCOORD')[:,0], Mesh.get('GCOORD')[:,1], Vel_x, Vel_y)
    plt.show()


    #plt.figure()
    #plt.triplot(Mesh.get('GCOORD')[:,0], Mesh.get('GCOORD')[:,1], Mesh.get('EL2NOD')[:,0:3])
    #plt.plot(Mesh.get('GCOORD')[:,0], Mesh.get('GCOORD')[:,1], 'o')
    #plt.show()
