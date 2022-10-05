Coupled diffusion problems
==========================
So far we have looked at problems that included a single unknown (e.g. temperature) at each node. Sometimes it is necessary to solve for multiple unknowns at each node simultaneously. This could just be multiple displacement or velocity components (as in mechanics), or this could be necessary when equations are strongly coupled.

One example, where a coupled system of equations needs to be solved are so-called Turing pattern :cite:`Turing1952`, which appear to emerge in many biological and chemical systems quasi spontaneously from near homogeneous starting conditions. A modern look at it can be found in the paper by :cite:`Maini2012`. Here we follow the implementation and problem description given in the FEM book by Guy Simpson :cite:`Simpson2017`. 

Governing equations
-------------------
We will explore the evolution of two chemical compounds, :math:`A` and :math:`B`, which affect each other. One is an 'activator' and one is an 'inhibitor'. Both compounds are also produced by a background rate and the their concentrations control each other. This can be expressed by a set of coupled equations.

.. math::
    :label: eq:turing_eqs
    
    \begin{align}
    \begin{split}
    \frac{\partial A}{\partial t} &= \nabla^2 A + \gamma (a - A +A^2 B)\\
    \frac{\partial B}{\partial t} &= d\nabla^2 B + \gamma (b - A^2 B)
    \end{split}
    \end{align}

In these equations, both compounds are produced at background rates :math:`\gamma a` and :math:`\gamma b`, :math:`A`decays with a rate :math:`\gamma A`, and reaction further produces the ' activator' :math:`A` resulting in a feedback that is described by the non-linear :math:`A^2 B` term, which couples the equations in a non-trivial way. 

FEM discretization
------------------
Using the same FD time discretization and using an implicit scheme, we can write these equations in FEM matrix form (just like we did for temperature in the previous example):

.. math::
    :label: eq:turing_eqs_fem

    \begin{align}
    \begin{split}
    (M(1+\Delta t\gamma) + \Delta tA_A)A^{n+1} = MA^n + \Delta tF_A^{n+1}\\
    (M + \Delta tA_B)B^{n+1} = MB^n + \Delta tF_B^{n+1}
    \end{split}
    \end{align}

with the matrices defined as:

.. math::
    :label: eq:turing_eqs_fem_matrix

    \begin{align}
    \begin{split}
    M &= \int_\Omega N_i N_j  d\Omega \ \ \ \ \ \ \ i,j=1,2,...,n\\
    A_A &= \int_\Omega \left ( \frac{\partial N_i}{\partial x}\frac{\partial N_j }{\partial x} + \frac{\partial N_i}{\partial y}\frac{\partial N_j }{\partial y} d\Omega \right ) d\Omega\ \ \ \ \ \ \ i,j=1,2,...,n\\
    A_B &= \int_\Omega d \left ( \frac{\partial N_i}{\partial x}\frac{\partial N_j }{\partial x} + \frac{\partial N_i}{\partial y}\frac{\partial N_j }{\partial y} d\Omega \right ) d\Omega\ \ \ \ \ \ \ i,j=1,2,...,n\\    
    F_A &= \int_\Omega \gamma N_i \left (a + (N_jA_j)^2N_jB_j      \right )  d\Omega \ \ \ \ \ \ \ i,j=1,2,...,n\\ 
    F_B &= \int_\Omega \gamma N_i \left (b - (N_jA_j)^2N_jB_j      \right )  d\Omega \ \ \ \ \ \ \ i,j=1,2,...,n\\ 
    \end{split}
    \end{align}

We can combine the matrices on the left-hand side and get the stifness matrices :math:`K_A` and :math:`K_B`.

.. math::
    :label: eq:turing_eqs_fem_matrix_2
    
    \begin{align}
    \begin{split}
    K_A A^{n+1} = M A^n + \Delta t F_A^{n+1}\\
    K_B B^{n+1} = M B^n + \Delta t F_B^{n+1}
    \end{split}
    \end{align}

Notice how we used a fully implicit scheme that evaluates the source terms at the new time level :math:`n+1`. We we will review this design decisions later.

One way of solving such a coupled problem is to have two unknowns per node (two so-called degrees of freedom). In our case, the complete element stiffness matrix would look like this:

.. math::
    :label: eq:turing_eqs_fem_matrix_3
    
    \begin{bmatrix}
    {K_A}_{11} & 0 & {K_A}_{12} & 0 & {K_A}_{13} & 0 \\
    0 & {K_B}_{11} & 0 & {K_B}_{12} & 0 & {K_B}_{13} \\
    {K_A}_{21} & 0 & {K_A}_{22} & 0 & {K_A}_{23} & 0 \\
    0 & {K_B}_{21} & 0 & {K_B}_{22} & 0 & {K_B}_{23} \\
    {K_A}_{31} & 0 & {K_A}_{32} & 0 & {K_A}_{33} & 0 \\
    0 & {K_B}_{31} & 0 & {K_B}_{32} & 0 & {K_B}_{33}
    \end{bmatrix}
    \begin{bmatrix}
    A_1^{n+1}\\
    B_1^{n+1}\\
    A_2^{n+1}\\
    B_2^{n+1}\\
    A_3^{n+1}\\
    B_3^{n+1}\\
    \end{bmatrix}
    = Rhs

The unknown concentrations of :math:`A` and :math:`B` or both showing up in the solution vector. What we also note is that the coupling is not that strong as there are no cross-terms. The equation for :math:`A_1^{n+1}`, the first row in :eq:`eq:turing_eqs_fem_matrix_3`, has zeros in the columns that operate on B. The coupling comes in through the source terms that appear on the Rhs (:eq:`eq:turing_eqs`).

Non-linear terms
----------------
Treating the non-linear source term in an implicit way requires iterations - because the term :math:`A^2B`cannot be directly incorporated into the stiffness matrix. We will use simple direct iterations, in which we solve for :math:`A` and :math:`B`, then update the source terms, and solve again until the solution is converged. Easiest to implement would be a non-linear iteration loop over the entire element loop but that would not be very smart as the stiffness matrix is not changing during iterations. Instead we first assembly the stiffness matrix and then have an iteration loop that updates the source term (right-hand side) until the solution is converged.


Python implementation
---------------------
We take our transient triangle solver from the previous excercise as a starting point and just modify the mesh and equation assembly, and add the non-linear iterations.

The model parameters are:
    * x,y have length 5
    * diffusivity for :math:`B` is :math:`d=20`
    * diffusivity for :math:`A` is :math:`1`
    * we use a triangle size of 0.005
    * and the other constants are:
        * :math:`\gamma=600`
        * :math:`a=0.05`
        * :math:`b=1`
        * time step is :math:`\Delta t=5e^{-4}`


Key pieces of the code are:

.. code-block:: python 
    :emphasize-lines: 10-13,30-31,34-35, 38,39,42,43
    :name: lst:2d-fem-turing_3
    :caption: Initial conditions
    
    nnodel = EL2NOD.shape[1]
    nel    = EL2NOD.shape[0]
    nnod   = GCOORD.shape[0]
    sdof   = nnod*2                 # two dof per node
    print(nnod, nel)

    # Initial conditions
    rng  = default_rng()
    vals = rng.standard_normal(nnod)
    A    = (a_coeff+b_coeff) + amp*vals
    vals = rng.standard_normal(nnod)
    B    = b_coeff/(a_coeff+b_coeff)**2 + amp*vals

This is how the initial conditions are set using random noise. Notice also how we define a variable :code:`sdof` that refers to the number of equations, while :code:`nnod` refers to the number of nodes in the mesh.

.. code-block:: python 
    :emphasize-lines: 10-13,30-31,34-35, 38,39,42,43
    :name: lst:2d-fem-turing_1
    :caption: Matrix assembly

    for t in range(0,nt):

    # Storage
    Rhs_all = np.zeros(sdof)
    I       = np.zeros((nel,2*nnodel*nnodel))
    J       = np.zeros((nel,2*nnodel*nnodel))
    K       = np.zeros((nel,2*nnodel*nnodel))
        for iel in range(0,nel):
            ECOORD  = np.take(GCOORD, EL2NOD[iel,:], axis=0 )
            Ael_A   = np.zeros((nnodel,nnodel))
            Ael_B   = np.zeros((nnodel,nnodel))
            RhsA_el = np.zeros(nnodel)
            RhsB_el = np.zeros(nnodel)
            
            for ip in range(0,nip):        
                # 1. update shape functions
                xi      = gauss[ip,0]
                eta     = gauss[ip,1]
                N, dNds = shapes_tri(xi, eta)
                
                # 2. set up Jacobian, inverse of Jacobian, and determinant
                Jac     = np.matmul(dNds,ECOORD) #[2,nnodel]*[nnodel,2]
                invJ    = np.linalg.inv(Jac)     
                detJ    = np.linalg.det(Jac)
                
                # 3. get global derivatives
                dNdx    = np.matmul(invJ, dNds) # [2,2]*[2,nnodel]
                
                # 4. compute element stiffness matrices
                Ael_A     = Ael_A + (np.outer(N,N)*(1+g_coeff*dt) +  dt*np.matmul(dNdx.T, dNdx))*detJ*weights[ip] 
                Ael_B     = Ael_B + (np.outer(N,N)        +  d_coeff*dt*np.matmul(dNdx.T, dNdx))*detJ*weights[ip] 
                
                # 5. assemble right-hand side
                RhsA_el     = RhsA_el + np.matmul(np.outer(N,N), np.take(A, EL2NOD[iel,:], axis=0 ))*detJ*weights[ip] 
                RhsB_el     = RhsB_el + np.matmul(np.outer(N,N), np.take(B, EL2NOD[iel,:], axis=0 ))*detJ*weights[ip] 
            
            # assemble coefficients
            I[iel,:]  =  np.concatenate((np.outer(2*EL2NOD[iel,:],np.ones(nnodel, dtype=int)).reshape(nnodel*nnodel),np.outer(2*EL2NOD[iel,:]+1,np.ones(nnodel, dtype=int)).reshape(nnodel*nnodel)))
            J[iel,:]  =  np.concatenate((np.outer(np.ones(nnodel, dtype=int),2*EL2NOD[iel,:]).reshape(nnodel*nnodel),np.outer(np.ones(nnodel, dtype=int),2*EL2NOD[iel,:]+1).reshape(nnodel*nnodel)))
            K[iel,:]  =  np.concatenate((Ael_A.reshape(nnodel*nnodel),Ael_B.reshape(nnodel*nnodel)))
            
            Rhs_all[2*EL2NOD[iel,:]]   += RhsA_el
            Rhs_all[2*EL2NOD[iel,:]+1] += RhsB_el

        A_all = csr_matrix((K.reshape(nel*2*nnodel*nnodel),(I.reshape(nel*2*nnodel*nnodel),J.reshape(nel*2*nnodel*nnodel))),shape=(sdof,sdof))


Notice how the node numbering is different to the equation numbering! Each node has two degrees of freedom, the concentration of :math:`A` and :math:`B`. All concentrations of :math:`A` are stored at :code:`2*EL2NOD[iel,:]` for the element, iel, and the concentrations of :math:`B` at :code:`2*EL2NOD[iel,:]+1`. This is a typical way of numbering the equations and we will see the same pattern when solving for viscous flow (when each node has two velocities).

.. code-block:: python 
    :emphasize-lines: 9-10,32,33,36,37, 40-44
    :name: lst:2d-fem-turing_2
    :caption: Non-linear iterations and solution

    # update right hand side in iterations

    iter = 0
    error = 10
    tol   = 0.001
    Conc_tmp = np.ones(sdof)*10
    iter_max = 20

    while error > tol:
        Tmp = Rhs_all.copy()
        iter += 1
        # loop over all elements and integrate Rhs
        for iel in range(0,nel):
            FA_el = np.zeros(nnodel)
            FB_el = np.zeros(nnodel)
            ECOORD  = np.take(GCOORD, EL2NOD[iel,:], axis=0 )

            for ip in range(0,nip):        
                # 1. update shape functions
                xi      = gauss[ip,0]
                eta     = gauss[ip,1]
                N, dNds = shapes_tri(xi, eta)
                
                # 2. set up Jacobian, inverse of Jacobian, and determinant
                Jac     = np.matmul(dNds,ECOORD) #[2,nnodel]*[nnodel,2]
                invJ    = np.linalg.inv(Jac)     
                detJ    = np.linalg.det(Jac)

                #3. integrate force vector
                Ai = np.dot(N,np.take(A, EL2NOD[iel,:], axis=0 ))
                Bi = np.dot(N,np.take(B, EL2NOD[iel,:], axis=0 ))
                FA_el     = FA_el + N*dt*g_coeff*(a_coeff+Bi*Ai**2)*detJ*weights[ip] # (dt*g_coeff*N*a_coeff+dt*g_coeff*N*np.dot(N,np.take(A, EL2NOD[iel,:], axis=0 ))**2*np.dot(N,np.take(B, EL2NOD[iel,:], axis=0 )))*detJ*weights[ip] 
                FB_el     = FB_el + N*dt*g_coeff*(b_coeff-Bi*Ai**2)*detJ*weights[ip] # (dt*g_coeff*N*b_coeff-dt*g_coeff*N*np.dot(N,np.take(A, EL2NOD[iel,:], axis=0 ))**2*np.dot(N,np.take(B, EL2NOD[iel,:], axis=0 )))*detJ*weights[ip] 
      
            # We don't have boundary conditions, as everything is zero flux      
            Tmp[2*EL2NOD[iel,:]]   += FA_el
            Tmp[2*EL2NOD[iel,:]+1] += FB_el
       
        # solve  system
        Conc  = spsolve(A_all,Tmp)        
        error = np.amax(np.absolute(Conc - Conc_tmp))/np.amax(np.absolute(Conc))
        Conc_tmp = Conc.copy()
        A     = Conc[0:sdof:2]
        B     = Conc[1:sdof:2]

        print(error, iter)
        if iter == iter_max:

            break

During the non-linear iterations, which happen inside the time loop, we keep changing the right-hand side and solve the system of equations. We do not re-assembly the stiffness matrix as it actually doesn't change during itreations. It actually also doesn't change during time steps and we could optimize the code further.


Excercise
---------

Try to puzzle a working code together and make sure that you understand the equation numbering and the way the iterations are done. If you get bored, stuck, or impatient: here is a link to the working code :download:`double_diff.py <python/double_diff.py>`.



.. only:: html

    Results of transient Turing problem.
 
    .. raw:: html
 
       <video width=100% autoplay muted controls loop>
       <source src="../_static/video/Turing.mp4" type="video/mp4">
          Your browser does not support HTML video.
       </video>





Is FEM a good idea for solving this problem?
---------------------------------------------
Ok, we got it to work but the solution seems super slow. We iterate to get a consistent right-hand side and solve the equations on an unstructured mesh although we are not resolving any complex geometries. 

A closer look at the equations also reveals that we are solving a problem that is dominated by the source terms (the reactions rates). Such problems are always a bit tricky to solve as the reactions often take please at faster rates than the diffusion process. As a consequence, we actually struggle to take big time steps and our fully implicit scheme is not efficient.

A much more effective way to solve the equations is to use explicit methods that are time-step limited (remember the FDM lectures?) but are also very efficient as they do not involve solving a system of equations. 

To illustrate this point, we look at a Turing problem that is described `here <https://blogs.mathworks.com/graphics/2015/03/16/how-the-tiger-got-its-stripes/>`_ and `here <http://www.karlsims.com/rd.html>`_ . Read those linked documents - they provide a really nice introduction to Turing pattern!

Alright, back to explicit FDM, let's see how badly FEM on unstructured meshes "loses" against FDM on a structured mesh. Have a look at the notebook and complete it! Try also other variations of "kill" and "feed" parameters.

.. toctree::
    :maxdepth: 2

    jupyter/turing_fdm.ipynb


