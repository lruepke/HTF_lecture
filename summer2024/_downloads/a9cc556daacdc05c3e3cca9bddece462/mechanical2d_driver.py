import numpy as np
import argparse
import os
import matplotlib.pyplot as plt
import triangle as tr
from mechanical2d import *
from write_triangle_files import *
import matplotlib.tri as tri

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
        segment_markers.extend(1001*np.ones(len(Tmp)))

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
    #    Node_ids    = np.hstack((Node_ids, 0*np.ones(new_nodes.shape[0])))
    #Node_ids    = np.hstack((Node_ids, np.zeros(new_nodes.shape[0], dtype=int)))
    Node_ids = np.hstack((Node_ids, np.zeros(new_nodes.shape[0]))).astype(int)
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
    D           = np.array([1e0,  1e3])      #Viscosity
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
    el_sizes    = (1e-4, 1e-4)
    el_ids      = (0,1)   
    write_files = True

    # and make it
    Mesh = make_mesh(x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, c_inc=c_inc, no_pts=no_pts,radius=radius, el_sizes=el_sizes, el_ids=el_ids)

    #write out number of nodes and elements
    print(f"Number of nodes: {Mesh.get('nnod')}")
    print(f"Number of elements: {Mesh.get('nel')}")

    if write_files:
        
        model_name = 'model'
        write_mesh_files(Mesh, model_name)

    #set boundary conditions
    #note that Bc_ind points to degree of freedeom and not node number!
    Bc_ind  = np.where(np.isin(Mesh.get('Node_ids'), [101, 102, 103, 104]))[0]
    Bc_val  = np.hstack((0.5*Mesh.get('GCOORD')[Bc_ind,0],  -0.5*Mesh.get('GCOORD')[Bc_ind,1]))
    Bc_ind  = np.hstack((2*Bc_ind, 2*Bc_ind+1))

    # and add to Mesh structure
    Mesh['Bc_ind'] = Bc_ind
    Mesh['Bc_val'] = Bc_val

    Materials = dict(D=D, Rho=Rho, G=G)

    Solver = dict(nip=6)

    # now call solver
    Vel, Pressure = mechanical2d(Mesh=Mesh, Materials=Materials, Solver=Solver  )

    Vel_x = Vel[0::2]
    Vel_y = Vel[1::2]

    #PLotting
    #we make a new connectivity, so that each pressure value is at an independent node,
    #this way we can see the discontinuous pressure field that varies lineraly in each element.
    GC_BIG = np.vstack((Mesh.get('GCOORD')[Mesh.get('EL2NOD')[:,0:3],0].ravel(), Mesh.get('GCOORD')[Mesh.get('EL2NOD')[:,0:3],1].ravel())).reshape(2,-1).T
    EN_BIG = np.arange(3*Mesh.get('nel')).reshape(-1,3)

    plt.figure()
    levels = np.linspace(Pressure.min(), Pressure.max(), num=100)
    contours = plt.tricontourf(GC_BIG[:,0], GC_BIG[:,1], EN_BIG, Pressure, levels=levels, cmap='jet')
    plt.triplot(Mesh.get('GCOORD')[:,0], Mesh.get('GCOORD')[:,1], Mesh.get('EL2NOD')[:,0:3], color='lightgrey', linewidth=0.1)
    #plt.quiver(Mesh.get('GCOORD')[:,0], Mesh.get('GCOORD')[:,1], Vel_x, Vel_y)

    plt.colorbar(contours, label='Pressure Values')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')
    plt.show()

 




