import matplotlib.pyplot as plt
import numpy as np
import triangle as tr

# arrays to fill in with input
vertices = []
segments = []
regions = []

# make a box with given dims and place given attribute at its center
def make_box(x, y, w, h, attribute):

    i = len(vertices)
    
    vertices.extend([[x,   y],
                     [x+w, y],
                     [x+w, y+h],
                     [x,   y+h]])

    segments.extend([(i+0, i+1),
                     (i+1, i+2),
                     (i+2, i+3),
                     (i+3, i+0)])
    
    regions.append([x+0.01*w, y+0.01*h, attribute,0.02])

def make_inclusion(center_x, center_y, radius, points_inc, attribute):
    
    theta        = np.linspace(0,2*np.pi,points_inc, endpoint=False)
    xx           = np.cos(theta)
    yy           = np.sin(theta)

    i = len(vertices)   

    vertices.extend(np.array([center_x + radius*xx,center_y + radius*yy]).T)

    Tmp = np.array([np.arange(i, i+points_inc), np.arange(i+1, i+points_inc+1)]).T
    Tmp[-1,1] = i

    segments.extend(Tmp)
   
    regions.append([center_x, center_y, attribute,0.01])

# generate input    
make_box(-1.0, -1.0, 2, 2, 1)
make_inclusion(0, 0, 0.2, 30, 100)

A = dict(vertices=vertices, segments=segments, regions=regions)
B = tr.triangulate(A, 'pq33Aa')
   
tr.compare(plt, A, B)
plt.show()