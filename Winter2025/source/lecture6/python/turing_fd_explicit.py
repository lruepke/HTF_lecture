import numpy as np
from scipy.sparse.dia import dia_matrix
from tabulate import tabulate
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix
from scipy.linalg import cho_factor, cho_solve
import matplotlib as mpl
import matplotlib.pyplot as plt
import triangle as tr
import meshio
from numpy.random import default_rng
import time
from scipy import ndimage, misc
from matplotlib import animation

# FEM utilities (refactored 2025) - for validation only
from validation import check_cfl_condition, validate_parameters


# periodic boundary condition laplacian FD scheme
def my_laplacian(in_M):
    out = -in_M + 0.25*(np.roll(in_M,1,axis=1) + np.roll(in_M,-1,axis=1) + np.roll(in_M,1,axis=0) + np.roll(in_M,-1,axis=0))

    return out

# 0. model parameters
da      =  1   # diffusion A
db      = .5  # diffusion B
f       = .055 # feed rate
k       = .062 # kill rate

# Validate parameters
params = {
    'da': da,
    'db': db,
    'f': f,
    'k': k
}
validate_parameters(params)

#Mesh and initial conditions
nx      = 128
A       = np.ones((nx,nx))
B       = np.zeros((nx,nx))
B[np.ix_(np.arange(50,61),np.arange(50,71))] = 1
B[np.ix_(np.arange(60,81),np.arange(70,81))] = 1
Anew    = A.copy()
Bnew    = B.copy()

x = np.linspace(0, 127, nx)
y = np.linspace(0, 127, nx)
xv, yv = np.meshgrid(x, y)
# time stepping
dt      = .25
tottime = 5000
t       = 0

# Calculate grid spacing
dx = 127.0 / (nx - 1)

# Check CFL condition (CRITICAL for explicit method!)
# Use the larger diffusivity for the check
D_max = max(da, db)
is_stable, cfl, cfl_max = check_cfl_condition(dt, dx, D_max, method='explicit', ndim=2)
print(f"CFL number: {cfl:.4f} (max for stability: {cfl_max:.4f})")
if is_stable:
    print("✓ CFL condition satisfied - simulation should be stable")
else:
    print("✗ WARNING: CFL condition violated - simulation may be unstable!")
    print(f"  Reduce dt to < {cfl_max * dx**2 / D_max:.4f} for stability")

# run model
while t<tottime:
    #Anew = A + (da*ndimage.laplace(A) - A*np.square(B) + f*(1-A))*dt
    #Bnew = B + (db*ndimage.laplace(B) + A*np.square(B) - (k+f)*B)*dt
    
    Anew = A + (da*my_laplacian(A) - A*np.square(B) + f*(1-A))*dt
    Bnew = B + (db*my_laplacian(B) + A*np.square(B) - (k+f)*B)*dt
    
    A = Anew
    B = Bnew
    t = t+dt

# plotting
fig = plt.figure()
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig.add_axes([left, bottom, width, height]) 

cp = plt.contourf(x,y, B, cmap='gist_yarg')
plt.colorbar(cp)
plt.show()
