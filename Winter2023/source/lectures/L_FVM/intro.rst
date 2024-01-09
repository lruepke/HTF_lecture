.. include:: /include.rst_

.. _L_FVM_Intro:

FVM and OpenFoam
=====================================================

In the previous session, we have used openFoam to resolve hydrothermal convection in oceanic crust. In this lecture, we will have a deeper look into how openFoam solves the governing equations and how the numerical schemes work. 

OpenFoam uses the Finite Volume Method (FVM), which is a popular numerical technique used in computational fluid dynamics (CFD) and other areas of computational physics to solve partial differential equations (PDEs) that describe physical phenomena. It's particularly well-suited for problems involving fluid flow, heat transfer, and associated physical processes. Here's an overview of how FVM works:

Basic Concept
---------------

**Domain Discretization** The computational domain (the area or volume where the physical problem occurs) is divided into small, discrete control volumes (cells). These control volumes cover the entire domain without overlapping.

**Integral Form of PDEs** FVM operates on the integral form of the governing PDEs, as opposed to the differential form used in methods like the Finite Difference Method (FDM). The equations are integrated over each control volume.

**Flux Calculations** The key aspect of FVM is the calculation of fluxes of the conserved quantities like mass, momentum, energy across the faces of each control volume. These fluxes are used to determine how these quantities change within each control volume over time.

**Boundary Conditions** The method also incorporates boundary conditions, which define how the field behaves at the boundaries of the computational domain.

The general workflow in FVM is the following: Represent the variables like velocity, pressure, temperature at discrete locations in each control volume - either at the cell centers or at the cell faces. The integral form of the governing PDEs is discretized for each control volume. This involves expressing the rate of change of a quantity in a control volume in terms of the fluxes through its faces. The fluxes through the control volume faces are approximated using values of the variables at the centers or faces of the cells. Various schemes like upwind and central difference can be used for this approximation. The discretized equations form a system of algebraic equations. These equations are solved iteratively to find the field values (like pressure, velocity, etc.) throughout the domain. For steady-state problems, iterations continue until the solution converges. For transient problems, the solution is marched forward in time, updating the variables at each time step. One key advantage of the FVM is that it ensures conservation of physical quantities locally (in each control volume) and globally (across the entire domain); it can further handle complex geometries and irregular meshes. In summary, FVM is a versatile and powerful tool for solving complex physical problems numerically. 

Openfoam uses what is called a cell-centered finite volume method. This means that the variables are stored at the center of the cells. The fluxes are calculated at the cell faces. The fluxes are then used to update the variables at the cell centers.

Example problem: heat diffusion
-----------------------------------

Let's look at an example problem that solves for temperature diffusion, using the laplacianFoam solver. 

Geometry and governing equations
--------------------------------
The governing equation is the heat diffusion equation, which is given by:

.. math::
    :label: eq:fvm_laplacian_con
    
    \rho c_p\frac{\partial T}{\partial t} - \nabla \cdot k \nabla T = 0 

.. math::
    :label: eq:fvm_laplacian_dif
    
    \frac{\partial T}{\partial t} - \nabla \cdot D \nabla T =0

:math:`k` is the thermal conductivity and :math:`D` is the thermal diffusivity :math:`\frac{k}{\rho c_p}` . We here assume that :math:`D` is a constant value. See also :ref:`theory_heat_diffusion`. 

Below is a simple 2D setup for a heat diffusion test case. We will use it to learn about the details of how openFOAM's FV method works.

.. figure:: /_figures/boundaryConditions_FVM_regularBox.*
   :align: center
   :name: fvm_geometry_bcs_regularBox
   :figwidth: 100%

   Model geometry, boundary conditions and initial condition.

.. admonition:: What we need to solve ?

   Temperature field :math:`T` of the model region at time :math:`t`, 
   therefore :math:`T` is the only unknown variable of our problem. It is a scalar field at the center of each cell.


Domain discretization
---------------------------------------------------

Mesh structure and topology
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can use a simple blockMesh to discretize our domain. You can download a prepared case here: :download:`Regular box case <cases/regularBox.zip>`.

.. tab:: Regular mesh

   .. figure:: /_figures/mesh_FVM_regularBox.*
      :align: center
      :name: fig:polyMesh_regularBox
      :figwidth: 100%

      2D regular mesh structure and topology information of OpenFOAM polyMesh (See :ref:`/lectures/L_FVM/cases/jupyter/VisualizeResults.ipynb#1.-Read-and-plot-mesh-information` in the notebook).

.. tab:: Unstructured mesh

   .. figure:: /_figures/mesh_FVM_unstructured.*
         :align: center
         :name: fig:polyMesh_unstructured
         :figwidth: 100%

         2D unstructured mesh topology information (:download:`Regular box case <cases/unstructured.zip>`).


The mesh comprises cell centers, cell faces, and cell vertices. The cell centers are the points where the temperature field is stored. The cell faces are the boundaries of the cells. The cell vertices are the points that define the cell faces. Notice how the internal faces (the ones that are not boundary conditions) are numbered continuously starting from 0.

In addition, a connectivity is needed to describe the topology of the mesh. This connectivity describes which cells are connected to each other through the faces. In openFoam a highly efficent data structure is used to store the connectivity information. For each face, the **owner** cell and the **neighbour** cell are stored. The normal vector of face always points from the owner cell to the neighbour cell. The numbering of the vertices that make up the face is again done using the right-hand rule and the face normal always points from the cell with the lower index to the cell with the larger index. Boundary faces  have only an **owner** cell.

You can find details on the mesh description in the `OpenFoam manual <https://cfd.direct/openfoam/user-guide/v8-mesh-description/>`_. Have a look at it and check the structure of the files in :code:`constant/polymesh`. 


.. tip::
   Check your understanding of the openFoam mesh structure and topology by creating the blockMesh mesh for the case we have downloaded and visualizing it using paraview. To create the labels, use the GenerateIds filter in paraview. 


Spatial discretization: The diffusion term
---------------------------------------------------

The equation discretization step is performed over each element of the computational domain to yield an algebraic relation that connects the value of a variable in an element to the values of the variable in the neighboring elements :cite:`moukalled2016finite`. 

.. admonition:: Goal

   Transform partial differential :eq:`eq:fvm_laplacian_dif` into a set of algebraic equations: :math:`\mathbf{A} T = \mathbf{b}`.

1. Integrate :eq:`eq:fvm_laplacian_dif` over element :math:`C` (the green cell in :numref:`fig:polyMesh_regularBox`) that enables recovering its integral balance form as,

.. math::
    :label: eq:fvm_volume_int
    
    \iiint\limits_{V_C} \frac{\partial T}{\partial t} dV - \iiint\limits_{V_C} \nabla \cdot D\nabla T dV = 0

**Note** that the transient term :math:`\iiint\limits_{V_C} \frac{\partial T}{\partial t} dV` will be processed in the later section. Let's set it to zero for the moment and assume steady-state.

2. Transform the volume integral of the heat flux into a surface integral by applying the divergence theorem, 

.. math::
    :label: eq:fvm_surface_int
    
    - \iiint\limits_{V_C} \nabla \cdot D\nabla T dV = -\oint\limits_{\partial V_C} (D\nabla T)\cdot d\vec{S} = 0

The :eq:`eq:fvm_surface_int` is actually a heat balance over cell :math:`C`. 
It is basically the integral form of the original partial differential equation and involves **no approximation**.

The integrant is the **diffusive heat flux**, 

.. math::
   :label: eq:fvm_flux_D

   \vec{J}^{T,D} \equiv -D\nabla T 

3. Transform the surface integral in :eq:`eq:fvm_surface_int` as a summation over the control volume faces (**still no approximation**),

.. math::
    :label: eq:fvm_surface_sum
    
    \sum\limits_{f\sim faces(V_C)} \iint\limits_{f}\vec{J}^{T,D}_f \cdot d\vec{S}_f =0

Here :math:`f` denotes the boundary face of cell :math:`V_C`.

.. .. math::
..     :label: eq:fvm_surface_sum_expand
    
..     b_C = \color{red}{(D\nabla T)_{f_{23}} \cdot \vec{S}_{f_{23}}} + (D\nabla T)_{f_2} \cdot \vec{S}_{f_2} + (D\nabla T)_{f_3} \cdot \vec{S}_{f_3}+ (D\nabla T)_{f_4} \cdot \vec{S}_{f_4}

Now, the key problem is how to calculate flux integral on a face.

4. Evaluating flux integral on the faces

It's worth noting that there are only two kinds of cells in the discretized domain, one is **internal cell**, the other is **boundary cell**.
All faces of a internal cell are internal faces. The remain cells are boundary cells which contain at least one boundary face.

.. tab:: Internal cell

   For a specific internal cell :math:`C`, e.g. cell 12, shown in :numref:`fig:polyMesh_regularBox`, 
   the :eq:`eq:fvm_surface_sum` can be expanded as,

   .. math::
      :label: eq:fvm_surface_sum_expand

      \color{red}{\iint_{f_{23}}\vec{J}^{T,D}_{f_{23}} \cdot \vec{S}_{f_{23}}} + \iint_{f_{24}}\vec{J}^{T,D}_{f_{24}} \cdot \vec{S}_{f_{24}} + \iint_{f_{21}}\vec{J}^{T,D}_{f_{21}} \cdot \vec{S}_{f_{21}}  + \iint_{f_{5}}\vec{J}^{T,D}_{f_{5}} \cdot \vec{S}_{f_{5}} = 0

   Note that the surface normals :math:`\vec{S}` always point out of the cell; for example, :math:`\vec{S}_{f_{23}}` would have a positive and :math:`\vec{S}_{f_{21}}` a negative sign.

   4.1 Considering face :math:`f_{23}`, calculate the first term on the right hand side of the :eq:`eq:fvm_surface_sum_expand` (we have to introduce the **first approximation** at this step). Using a Gaussian quadrature the integral at the face :math:`f_{23}`, for example, becomes,

   .. math::
      :label: fvm_flux_Gaussian_integral

      \color{red}{\iint_{f_{23}}\vec{J}^{T,D}_{f_{23}} \cdot \vec{S}_{f_{23}}} = \iint_{f_{23}} (\vec{J}^{T,D}_{f_{23}} \cdot \vec{n}_{f_{23}}) dS_{f_{23}} \approx \color{orange}{\sum\limits_{ip\sim ip(f_{23})} (\vec{J}^{T,D} \cdot \vec{n}_{ip})\omega_{ip} S_{f_{23}}}

   where :math:`S_f` is the area of face :math:`f_{23}`, :math:`ip` refers to a integration point and :math:`ip(f_{23})` the number of integration points along surface :math:`f_{23}`, :math:`\omega_{ip}` is the integral weights. 

   .. figure:: /_figures/Fig5.2_RedBook.*
      :align: center
      :width: 100 %
      :name: fig:Fig5.2
   
      Surface integration of fluxes using (a) one integration point, (b) two integration points, and (c) three integration points :cite:`moukalled2016finite`. The total flux :math:`FluxT_f` is computed as linear combinations of the two cell values plus a source term.
   
   .. important:: 

      Only the **one integration point scheme**, i.e. :numref:`fig:Fig5.2` (a), is implemented in OpenFOAM (at least before version 8).
      Therefore the first keyword of Laplacian scheme in :code:`system/fvScheme` dictionary file of a case has only one option, it is :code:`Gauss`.

   4.2 Choose integral scheme or integral points

   To simply explain the calculation process/logic and to also be consistent with OpenFOAM, we here adopt an **one integration point** scheme with weight :math:`\omega = 1`, thus :eq:`fvm_flux_Gaussian_integral` becomes,

   .. math::
      :label: eq:fvm_surface_term_f1
      
      \color{orange}{\sum\limits_{ip\sim ip(f_{23})} (\vec{J}^{T,D} \cdot \vec{n}_{ip})\omega S_f } = -\left(D \frac{\partial T}{\partial x} \vec{i} + D\frac{\partial T}{\partial y} \vec{j} \right)_{f_{23}} \cdot \Delta y_{f_{23}} \vec{i} = -\color{blue}{\left( D\frac{\partial T}{\partial x} \right)_{f_{23}} \Delta y_{f_{23}}}

   where 

   * :math:`S_{f_{23}} = \Delta y_{f_{23}}`  is the area of face  :math:`f_{23}` by assuming :math:`\Delta z = 1`. 
   * :math:`\vec{S}_{f_{23}} = S_{f_{23}} \vec{n}_{f_{23}}` is the surface vector of face :math:`f_{23}`.
   * :math:`\vec{n}_{f_{23}} = \vec{i}` is the normal vector of the face :math:`f_{23}` directed out of the cell :math:`C` (see :numref:`fig:polyMesh_regularBox`).

   4.3 Calculate gradient of :math:`T` at face centroid, here introduce **the second approximation**, e.g. assuming linear variation of T and then the gradient term in the :eq:`eq:fvm_surface_term_f1` can be approximated as,

   .. math::

      -\color{blue}{\left( D\frac{\partial T}{\partial x} \right)_{f_{23}} \Delta y_{f_{23}}} = -D\frac{T_{13} - T_{C}}{x_{13} - x_C}\Delta y_{f_{23}} = -\frac{D \Delta y_{f_{23}}}   {\delta x_{f_{23}}}   (T_{13} - T_{C})

   where :math:`\delta x_{f_{23}}` represents the distance between center of cells who share face :math:`f_{23}`, they are cell 12 and cell 13.

   Now let's look at the western face :math:`f_{21}` . Using a single integration point, we can again write:
 
   .. math::
      :label: eq:fvm_surface_term_f1_west
      
      \color{orange}{\sum\limits_{ip\sim ip(f_{21})} (\vec{J}^{T,D} \cdot \vec{n}_{ip})\omega S_f } = -\left(D \frac{\partial T}{\partial x} \vec{i} + D\frac{\partial T}{\partial y} \vec{j} \right)_{f_{21}} \cdot \Delta y_{f_{21}} \cdot-\vec{i} = \color{blue}{\left( D\frac{\partial T}{\partial x} \right)_{f_{21}} \Delta y_{f_{21}}}

   Notice how the surface normal now has the opposite sign!

   We continue to spell out the flux across face :math:`f_{21}` as:

    .. math::

      \color{blue}{\left( D\frac{\partial T}{\partial x} \right)_{f_{21}} \Delta y_{f_{21}}} = D\frac{T_{C} - T_{11}}{x_{11} - x_C}\Delta y_{f_{21}} = -\frac{D \Delta y_{f_{21}}}   {\delta x_{f_{21}}}   (T_{11} - T_{C})


   Do the same thing for the remain faces (:math:`f_5, f_{24}`).
   


   4.3 Construct coefficients of a specific cell :math:`C_{12}` 

   Now we get all the coefficients, thus the :eq:`eq:fvm_surface_sum_expand` can be discretized and expressed as a matrix form,

   .. math::
      :label: eq:fvm_matrix_form_internalCell
      
      a_{C} T_C + a_{F_{23}} T_{F_{23}} + a_{F_{24}} T_{F_{24}} + a_{F_{21}} T_{F_{21}} + a_{F_{5}} T_{F_{5}} = 0


   with

   .. math::

      \begin{eqnarray}   
      a_{23} & = &  -\frac{D \Delta y_{f_{23}}}   {\delta x_{f_{23}}} \\
      a_{24} & = &  -\frac{D \Delta y_{f_{24}}}   {\delta x_{f_{24}}} \\
      a_{21} & = &  -\frac{D \Delta y_{f_{21}}}   {\delta x_{f_{21}}} \\
      a_{5} & = &  -\frac{D \Delta y_{f_{5}}}   {\delta x_{f_{5}}} \\
      a_{C} & = & -(a_{23} + a_{24} + a_{21}  + a_{5}) \\
      \end{eqnarray}
   
   .. math:: 
      :label: eq:fvm_matrix_form_internalCell12

      \mathbf{A}_{N\times N}[12,:] \mathbf{T}_{N\times 1} =  0
   
.. tab:: Boundary cell

   For a specific boundary cell, e.g. **cell 19** shown in :numref:`fig:polyMesh_regularBox`, the :eq:`eq:fvm_surface_sum` can be rewritten as,

   .. math::
      :label: eq:fvm_surface_sum_expand_boundary

      \iint_{f_{18}}\vec{J}^{T,D}_{f_{18}} \cdot \vec{S}_{f_{18}} + \iint_{f_{35}}\vec{J}^{T,D}_{f_{35}} \cdot \vec{S}_{f_{35}} +\iint_{f_{37}}\vec{J}^{T,D}_{f_{37}} \cdot \vec{S}_{f_{37}} + \color{red}{\iint_{f_{51}}\vec{J}^{T,D}_{f_{51}} \cdot \vec{S}_{f_{51}}} = 0
   
   The first three terms are normal fluxes across internal faces, so there is nothing special and we can just calculate them following steps for internal faces as before. The special thing is the boundary face :math:`f_{52}`, we can call it :math:`f_b` for a general case. **Now the problem is how to evaluate flux on the boundary face** :math:`f_b`.

   The fluxes on the interior faces are discretized as before, while **the boundary flux is discretized with the aim of constructing a linearization with respect to** the cell field :math:`T_C`, e.g. :math:`T_{C_{19}}` of cell :math:`C_{19}` shown in :numref:`fig:polyMesh_regularBox`, thus 

   .. math::
      :label: eq:fvm_surface_flux_boundary

      \color{red}{\iint_{f_{b}}\vec{J}^{T,D}_{f_{b}} \cdot d\vec{S}_{f_{b}}} = a_{F_b} T_C + c_{F_b}

   All right, **now out goal is to determine the coefficients of** :math:`a_{F_b}` and :math:`c_{F_b}` according to the boundary conditions. Generally there are two basic kinds of boundary conditions, they are **Dirichlet boundary condition** (also called the first type or fixed value boundary condition) and **Von Neumann boundary condition** (all called the second type or fixed flux boundary condition), respectively.

   4.1 Fixed value boundary condition

   Fixed value means the value of field :math:`T` is given on the boundary face is give, i.e. :math:`\color{purple}{T_{f_b}}` is given. There fore we can calculate flux on boundary face, :math:`f_{51}` for example,

   .. math::
      :label: eq:fvm_flux_boundary

      \color{red}{\iint_{f_{51}}\vec{J}^{T,D}_{f_{51}} \cdot \vec{S}_{f_{51}}} = \iint_{f_{51}} (\vec{J}^{T,D}_{f_{51}} \cdot \vec{n}_{f_{51}}) dS_{f_{51}} \approx \color{orange}{\sum\limits_{ip\sim ip(f_{51})} (\vec{J}^{T,D} \cdot \vec{n}_{ip})\omega_{ip} S_{f_{51}}}

   Here we still use one integration point scheme to explain the calculation process,

   .. math::
      :label: eq:fvm_flux_Gaussian_integral_boundary

      \color{orange}{\sum\limits_{ip\sim ip(f_{51})} (\vec{J}^{T,D} \cdot \vec{n}_{ip})\omega S_f } = -\left(D \frac{\partial T}{\partial x} \vec{i} + D\frac{\partial T}{\partial y} \vec{j} \right)_{f_{51}} \cdot \Delta y_{f_{51}} \vec{i} = -\color{blue}{\left( D\frac{\partial T}{\partial x} \right)_{f_{51}} \Delta y_{f_{51}}}

   Now let's calculate the gradient on the boundary face,

   .. math::
      :label: eq:fvm_gradient_boundary

      -\color{blue}{\left( D\frac{\partial T}{\partial x} \right)_{f_{51}} \Delta y_{f_{51}}} = -D\frac{\color{purple}{T_{f_{51}}} - T_{C}}{x_{f_{51}} - x_C}\Delta y_{f_{51}} = -\frac{D \Delta y_{f_{51}}}{\delta x_{f_{51}}}(\color{purple}{T_{f_{51}}} - T_{C}) = a_{f_{51}}(-T_{C}) +  c_{f_{51}}

   Substituting :eq:`eq:fvm_gradient_boundary` back to  :eq:`eq:fvm_flux_Gaussian_integral_boundary`, :eq:`eq:fvm_flux_boundary` and  :eq:`eq:fvm_surface_flux_boundary`, we can get Laplacian discretization coefficients :math:`a_{f_{51}}` and :math:`c_{f_{51}}`,

   .. math::
      :label: eq:fvm_laplacian_coeff_boundary_fixedvalue

      \begin{eqnarray}
      a_{F_{51}} &=& -D\frac{\Delta y_{f_{51}}}{\delta x_{f_{51}}}\\
      c_{F_{51}} &=& a_{F_{51}}\color{purple}{T_{f_{51}}}
      \end{eqnarray}

   Notice how :math:`\color{purple}{T_{f_{51}}}` is known (a boundary condition), so that the term :math:`-a_{F_{51}}\color{purple}{T_{f_{51}}}` will end up on the right-hand side and not in the coefficient matrix.

   4.2 Fixed flux boundary condition

   **Fixed flux** means **normal gradient to the boundary face** of field :math:`T` is given, i.e. :math:`\color{purple}{\vec{J}^{T,D}_{f_{b}} \cdot \vec{n}_{f_{b}} = q_{f_b}}` is given. There fore we can calculate flux on boundary face, :math:`f_{51}` for example,

   .. math::
      :label: eq:fvm_fixedFlux_boundary

      \color{red}{\iint_{f_{51}}\vec{J}^{T,D}_{f_{51}} \cdot \vec{S}_{f_{51}}} = \iint_{f_{51}} (\vec{J}^{T,D}_{f_{51}} \cdot \vec{n}_{f_{51}}) dS_{f_{51}} = \iint_{f_{51}} \color{purple}{q_{f_{51}}} dS_{f_{51}} = \color{purple}{q_{f_{51}}} S_{f_{51}} = \color{purple}{q_{f_{51}}} \Delta y_{f_{51}}

   Substituting :eq:`eq:fvm_fixedFlux_boundary` back to :eq:`eq:fvm_surface_flux_boundary`, we can get Laplacian discretization coefficients :math:`a_{f_{51}}` and :math:`c_{f_{51}}`,

   .. math::
      :label: eq:fvm_laplacian_coeff_boundary_fixedflux

      \begin{eqnarray}
      a_{F_{51}} &=& 0 \\
      c_{F_{51}} &=& \color{purple}{q_{f_{51}}} \Delta y_{f_{51}}
      \end{eqnarray}

   4.3 Construct coefficients of the boundary cell :math:`C_{19}`

   Do the same thing as 4.1 or 4.2 (depends on the type of boundary condition) for all boundary faces of a boundary cell :math:`C` and calculate the coefficients of :math:`a_{F_b}` and :math:`c_{F_b}`, then the :eq:`eq:fvm_surface_sum_expand_boundary` can be expressed as,

   .. math::
      :label: eq:fvm_matrix_form_boundaryCell
      
        \begin{eqnarray}
      b_{C_{19}} & =  & a_{F_{18}} (T_{F_{18}} - T_C) + a_{F_{35}} (T_{F_{35}} - T_C) + a_{F_{37}} (T_{F_{37}} - T_C) + a_{F_{51}} (-T_C) +c_{F_{51}} \\
         & = & \sum\limits_{f\in[18, 35, 37]} a_{F_f} T_{F_f} + \left(-\sum\limits_{f\in[18, 35, 37, 51]} a_{F_f} \right) T_C + \sum\limits_{f\in[51]} c_{F_f} \\
         & = & \sum\limits_{f\sim nb(\text{internal face})} a_{F_f} T_{F_f} + \left(-\sum\limits_{f\sim nb(\text{all face})} a_{F_f} \right) T_C + \sum\limits_{f\sim nb(\text{boundary face})} c_{F_f}
      \end{eqnarray}

   .. math:: 
      :label: eq:fvm_matrix_form_boundaryCell19

      \mathbf{A}_{N\times N}[19,:] \mathbf{T}_{N\times 1} =  B_{19}

   where :math:`B_{19} = b_{C_{19}} - \sum\limits_{f\sim nb(\text{boundary face})} c_{F_f}`.

5. Construct general matrix form of the :eq:`eq:fvm_surface_sum`

.. math::
   :label: eq:fvm_matrix_form0

   \begin{eqnarray}
   -\sum\limits_{f\sim nb(C)} (D\nabla T)_{F_f}\cdot \vec{S}_{F_f} & = &  \sum\limits_{f\sim \text{internal faces}(C)} a_{F_f} T_{C} + \left(-\sum\limits_{f\sim nb(C)} a_{F_f} T_{F_f} \right) + \sum\limits_{f\sim \text{boundary faces}(C)} c_{F_{F_f}} \\
   & = & a_CT_C - \sum\limits_{f\sim \text{internal faces}(C)} a_{F_f} T_{F_f} + \sum\limits_{f\sim \text{boundary faces}(C)} c_{F_{f}}
   \end{eqnarray}

Matrix form of the :eq:`eq:fvm_surface_sum` can be written as,

.. math::
   :label: fvm_matrix_form

   \mathbf{A}_{N \times N} \mathbf{T}_{N\times1} = \mathbf{B}_{N\times1}

here :math:`N` is the number of cell, :math:`\mathbf{A}` is a sparse matrix of coefficients, 
:math:`\mathbf{T}` is cell-centered temperature field.
The matrix is visualized as :numref:`fig:fvm_matrix`.
It should be noted that **the coefficients matrix of Laplacian term is a symmetric matrix and the boundary conditions only affect diagonal entry of the coefficients matrix**. Further, only the fixed value boundary condition affects the matrix and the fixed flux boundary condition affects the RHS(Right hand of side) :math:`B`.

.. tab:: Regular mesh 

   .. figure:: /_figures/matrix_FVM_regularBox.*
      :align: center
      :name: fig:fvm_matrix
      :figwidth: 100%

      Visualization of :math:`\mathbf{A}_{N \times N} \mathbf{T}_{N\times1} = \mathbf{b}_{N\times1}`. 
      The coefficients of the selected cell :math:`C` (:numref:`fig:polyMesh_regularBox`) are marked by green rectangle.

.. tab:: Unstructured mesh 

   .. figure:: /_figures/matrix_FVM_unstructured.*
      :align: center
      :name: fig:fvm_matrix_unstructured
      :figwidth: 100%

      Visualization of :math:`\mathbf{A}_{N \times N} \mathbf{T}_{N\times1} = \mathbf{b}_{N\times1}`. 
      The coefficients of the selected cell :math:`C` (:numref:`fig:polyMesh_regularBox`) are marked by green rectangle.


Temporal Discretization: The Transient Term
------------------------------------------------------

After spatial discretization, the :eq:`eq:fvm_volume_int` can be expressed as,

.. math::
   :label: eq:fvm_diffusionEq_spatial_dis

   \frac{\partial T}{\partial t} V_C - L(T_C^t) =0

where :math:`V_C` is the volume of the discretization cell and :math:`L(T_C^t)` is the spatial discretization operator expressed at some reference time :math:`t`, which can be written as algebraic form (see also :eq:`eq:fvm_matrix_form0`),

.. math::
   :label: eq:fvm_spatial_dis_algebraic

   L(T_C^t) = a_C T_C^t - \sum\limits_{F\sim NB(C)} a_F T_F^t + \sum\limits_{F\sim NB(C)} c_F

where :math:`a_C` is the diagonal coefficients of the matrix, :math:`a_F` is the off-diagonal coefficients, and the :math:`c_F` is the source coefficients as right hand side of matrix of system. 

.. Tip::

   For specific cell :math:`C`, 
   * :math:`a_F` has only contributed from internal faces. 
   * :math:`a_C` is the negative summation of :math:`a_F`, contributed from all faces, but equation of the :math:`a_F` for a boundary faces (:eq:`eq:fvm_laplacian_coeff_boundary_fixedvalue` and :eq:`eq:fvm_laplacian_coeff_boundary_fixedflux`) is a little bit different from internal face. 
   
   * :math:`c_F` only comes from boundary faces of cell :math:`C` if it has boundary face, see also :eq:`eq:fvm_laplacian_coeff_boundary_fixedvalue` and :eq:`eq:fvm_laplacian_coeff_boundary_fixedflux`. :math:`c_F` will contribute to RHS (:math:`\mathbf{B}`) of the algebraic :eq:`fvm_matrix_form`.
   
Let's do the temporal discretization for the first term of :eq:`eq:fvm_diffusionEq_spatial_dis`,

.. math:: 

   \frac{T^{t+\Delta t/2} - T^{t - \Delta t/2}}{\Delta t} V_C + L(T_C^t) = 0

To derive the full discretized equation, an interpolation profile expressing the face values at (:math:`t-\Delta t/2`) and (:math:`t+\Delta t`) in terms of the element values at (:math:`t`), (:math:`t-\Delta t`), etc., is needed.Independent of the profile used, the flux will be linearized based on old and new values as,

.. math::
   :label: eq:fvm_temporal_linear 

   b_C = FluxC T_C + FluxC^o T_C^o + FluxV

where the superscript :math:`^{o}` refers to old values. With the format of the linearization, 
the coefficient :math:`FluxC` will be **added to diagonal element** and the coefficient :math:`FluxC^o T_C^o + FluxV` will be **added to the source or RHS** (:math:`B`).

First order implicit Euler scheme
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::
   :label: eq:fvm_firstOrder_Euler_imp

   \frac{T^t - T^{t-\Delta t}}{\Delta t} V_C - L(T^{t}) = 0

Therefore the coefficients will be 

.. math::
   :label: eq:fvm_firstOrder_Euler_imp_coeff

   FluxC = \frac{V_C}{\Delta t}, ~ FluxC^o = - \frac{V_C}{\Delta t}, ~ FluxV =0

First order explicit Euler scheme
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math:: 
   :label: eq:fvm_firstOrder_Euler_exp

   \frac{T^{t+\Delta t} - T^t}{\Delta t} V_C - L(T^{t}) = 0

Note that now the new time is at :math:`t+\Delta t` and that the spatial operator (:math:`L`) of :eq:`eq:fvm_firstOrder_Euler_exp` has to be evaluated at time :math:`t`. 
Therefore, in order to find the value of :math:`T` at time :math:`t+\Delta t`, we don't need to solve a set of linear algebraic equations,

.. math::
   :label: eq:fvm_firstOrder_Euler_exp_coef

   T^{t+\Delta t} = \frac{\Delta t}{V_C} L(T^{t}) + T^t


Building the matrix: Mesh connectivity and internal data storage
----------------------------------------------------------------

Remember the mesh structure discussed above. The mesh connectivity is stored in the :code:`polyMesh` folder, which containes the :code:`neighbour` and :code:`owner` files. The normal vector of the face always points from the owner cell to the neighbour cell. The numbering of the vertices that make up the face is again done using the right-hand rule and the face normal always points from the cell with the lower index to the cell with the larger index. Boundary faces have only an owner cell.

Is is straightforward to find the indices for setting up the flux balance for each cell, as we did in :eq:`eq:fvm_matrix_form_internalCell`?

No, because the owner/neighbor data structure does not provide the indices of all the cells that are connected to a given cell. Rather the owner/neighbor data structure is "face-centered" in that for each face, we know which cells are on each side. Hence, we can easily loop over all faces and compute the fluxes but we cannot easily loop over cells and compute the fluxes to/from the neighboring cells. 

.. tip::
   Check out how openfoam makes use of this face-centered logic to compute the average gradient of a field in the :code:`gaussGrad` function in :code:`src/finiteVolume/finiteVolume/gradSchemes/gaussGrad/gaussGrad.C` (see below). The forAll loop is on the faces and the owner and neighbour functions are used to get the indices of the cells that are connected to a given face.

   .. code-block:: foam 
      :caption: src/finiteVolume/finiteVolume/gradSchemes/gaussGrad/gaussGrad.C
      :emphasize-lines: 82-88
      :linenos:

      /*---------------------------------------------------------------------------*\
        =========                 |
        \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
         \\    /   O peration     | Website:  https://openfoam.org
          \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
           \\/     M anipulation  |
      -------------------------------------------------------------------------------
      License
      This file is part of OpenFOAM.

      OpenFOAM is free software: you can redistribute it and/or modify it
      under the terms of the GNU General Public License as published by
      the Free Software Foundation, either version 3 of the License, or
      (at your option) any later version.

      OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
      ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
      FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
      for more details.

      You should have received a copy of the GNU General Public License
      along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

      \*---------------------------------------------------------------------------*/

      #include "gaussGrad.H"
      #include "extrapolatedCalculatedFvPatchField.H"

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

      template<class Type>
      Foam::tmp
      <
         Foam::GeometricField
         <
            typename Foam::outerProduct<Foam::vector, Type>::type,
            Foam::fvPatchField,
            Foam::volMesh
         >
      >
      Foam::fv::gaussGrad<Type>::gradf
      (
         const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
         const word& name
      )
      {
         typedef typename outerProduct<vector, Type>::type GradType;

         const fvMesh& mesh = ssf.mesh();

         tmp<GeometricField<GradType, fvPatchField, volMesh>> tgGrad
         (
            new GeometricField<GradType, fvPatchField, volMesh>
            (
                  IOobject
                  (
                     name,
                     ssf.instance(),
                     mesh,
                     IOobject::NO_READ,
                     IOobject::NO_WRITE
                  ),
                  mesh,
                  dimensioned<GradType>
                  (
                     "0",
                     ssf.dimensions()/dimLength,
                     Zero
                  ),
                  extrapolatedCalculatedFvPatchField<GradType>::typeName
            )
         );
         GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad.ref();

         const labelUList& owner = mesh.owner();
         const labelUList& neighbour = mesh.neighbour();
         const vectorField& Sf = mesh.Sf();

         Field<GradType>& igGrad = gGrad;
         const Field<Type>& issf = ssf;

         forAll(owner, facei)
         {
            GradType Sfssf = Sf[facei]*issf[facei];

            igGrad[owner[facei]] += Sfssf;
            igGrad[neighbour[facei]] -= Sfssf;
         }

         forAll(mesh.boundary(), patchi)
         {
            const labelUList& pFaceCells =
                  mesh.boundary()[patchi].faceCells();

            const vectorField& pSf = mesh.Sf().boundaryField()[patchi];

            const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

            forAll(mesh.boundary()[patchi], facei)
            {
                  igGrad[pFaceCells[facei]] += pSf[facei]*pssf[facei];
            }
         }

         igGrad /= mesh.V();

         gGrad.correctBoundaryConditions();

         return tgGrad;
      }


      template<class Type>
      Foam::tmp
      <
         Foam::GeometricField
         <
            typename Foam::outerProduct<Foam::vector, Type>::type,
            Foam::fvPatchField,
            Foam::volMesh
         >
      >
      Foam::fv::gaussGrad<Type>::calcGrad
      (
         const GeometricField<Type, fvPatchField, volMesh>& vsf,
         const word& name
      ) const
      {
         typedef typename outerProduct<vector, Type>::type GradType;

         tmp<GeometricField<GradType, fvPatchField, volMesh>> tgGrad
         (
            gradf(tinterpScheme_().interpolate(vsf), name)
         );
         GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad.ref();

         correctBoundaryConditions(vsf, gGrad);

         return tgGrad;
      }


      template<class Type>
      void Foam::fv::gaussGrad<Type>::correctBoundaryConditions
      (
         const GeometricField<Type, fvPatchField, volMesh>& vsf,
         GeometricField
         <
            typename outerProduct<vector, Type>::type, fvPatchField, volMesh
         >& gGrad
      )
      {
         typename GeometricField
         <
            typename outerProduct<vector, Type>::type, fvPatchField, volMesh
         >::Boundary& gGradbf = gGrad.boundaryFieldRef();

         forAll(vsf.boundaryField(), patchi)
         {
            if (!vsf.boundaryField()[patchi].coupled())
            {
                  const vectorField n
                  (
                     vsf.mesh().Sf().boundaryField()[patchi]
                  / vsf.mesh().magSf().boundaryField()[patchi]
                  );

                  gGradbf[patchi] += n *
                  (
                     vsf.boundaryField()[patchi].snGrad()
                  - (n & gGradbf[patchi])
                  );
            }
         }
      }


      // ************************************************************************* //




OpenFoam uses a matrix format that is tightly integrated with the "face-centered" and the owner/neighbor structure. It is called the LDU matrix format. This format is a way of storing sparse matrices. A sparse matrix is one in which most elements are zero, which is a common case in CFD (check the example above). The LDU format is specifically tailored to store and manipulate these types of matrices efficiently.

**Components of LDU:**
   - L (Lower triangle): This part of the matrix stores the coefficients that are below the main diagonal.
   - D (Diagonal): This represents the main diagonal of the matrix. In many algorithms, the diagonal elements play a crucial role and are treated separately for reasons of numerical stability and efficiency.
   - U (Upper triangle): This stores the coefficients above the main diagonal.

The LDU format is memory-efficient because it only stores non-zero elements. This is particularly beneficial in CFD where matrices can be very large, but only a small fraction of the elements are non-zero. OpenFOAM uses this structure to perform matrix operations like multiplication, addition, and especially the solving of linear systems, which is crucial in the iterative methods used in CFD simulations.

The LDU format uses scalar fields representing the diagonal, upper, and lower coefficients. The diagonal coefficients are indexed using the cell index. The upper and lower coefficients use face indices. In order to get the indices to match, mapping arrays between the face and cell indices are provided.  The :code:`lowerAddr()`, and :code:`upperAddr()` functions of the lduAdressing class provide the owner and neighbor cell indices for each face. 

A nice summary of the LDU format can be found in the `OpenFoam Wiki <https://openfoamwiki.net/index.php/OpenFOAM_guide/Matrices_in_OpenFOAM>`_.

.. admonition:: Summary

   In summary, the LDU format is a way of storing sparse matrices. It is specifically tailored to store and manipulate these types of matrices efficiently. OpenFOAM uses this structure to perform matrix operations like multiplication, addition, and especially the solving of linear systems, which is crucial in the iterative methods used in CFD simulations. Most of this happens under the hood but if you ever want to perform an operation on all rows of an lduMatrix, you need to know abou lduAdressing.






