.. include:: /include.rst_

.. _L_FVM_Exercise1:

.. warning::

   Before we start the exercise 1, please download (:download:`test_laplacianFoam <cases/test_laplacianFoam.zip>`) the modified Laplacian solver first.
   Then put it in the shared folder and compile (:code:`wmake`) it in the docker container.
   This modified solver writes out detailed information on matrix coefficients into the log file, so that it is possible to check the derived values against the computed values. 

==============
Exercise 1
==============


Theory
===========


.. admonition:: Goals

   - explore FVM theory
   - link FVM theory and OpenFOAM solver implementation
   - step by step exploring details of the solver and solving process


Step 1, Geometric and physical modeling
---------------------------------------------------

.. math::
    :label: eq:fvm_laplacian_con
    
    \rho c_p\frac{\partial T}{\partial t} - \nabla \cdot k \nabla T = 0 

.. math::
    :label: eq:fvm_laplacian_dif
    
    \frac{\partial T}{\partial t} - \nabla \cdot D \nabla T =0

:math:`k` is the thermal conductivity and :math:`D` is the thermal diffusivity :math:`\frac{k}{\rho c_p}` . We here assume that :math:`D` is a constant value.
see also :ref:`theory_heat_diffusion`. 

Below is a simple 2D setup for a heat diffusion test case. We will use it to learn about the details of how openFOAM's FV method works.

.. figure:: /_figures/boundaryConditions_FVM_regularBox.*
   :align: center
   :name: fvm_geometry_bcs_regularBox
   :figwidth: 100%

   Model geometry, boundary conditions and initial condition.

.. admonition:: What we need to solve ?

   Temperature field :math:`T` of the model region at time :math:`t`, 
   therefore :math:`T` is the only unknown variable of our problem.


Step 2, Domain discretization
---------------------------------------------------

Mesh structure and topology
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a 2D regular box mesh using :code:`blockMesh` utility (:download:`Regular box case <cases/regularBox.zip>`).
Therefore, the mesh topology shown below is the same as OpenFOAM polyMesh.

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

.. tip::

   All internal faces have and only have two attributes, **owner** and **neighbour**, 
   denote the cell indices who share the face. 
   While the boundary faces only have **owner** attribute.
   The normal vector of face always points from the **owner** cell to the **neighbour** cell. The numbering of the vertices that make up the face is again done using the right-hand rule and the face normal always points from the cell with the lower index to the cell with the larger index.

   You can find details on the mesh description in the `OpenFoam manual <https://cfd.direct/openfoam/user-guide/v8-mesh-description/>`_. Have a look at it and check the structure of the files in :code:`constant\polymesh`. 


.. tip::
   Check your understanding of the openFoam mesh structure and topology by creating the blockMesh mesh for the case we have downloaded and visualizing it using paraview. To create the labels, use the GenerateIds filter in paraview. 


Step 3, Spatial discretization: The diffusion term
---------------------------------------------------

The equation discretization step is performed over each element of the computational domain to yield an algebraic relation that connects the value of a variable in an element to the values of the variable in the neighboring elements :cite:`moukalled2016finite`. 

.. admonition:: Goal

   Transform partial differential :eq:`eq:fvm_laplacian_dif` into a set of algebraic equations: :math:`\mathbf{A}[T] = \mathbf{b}`.

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


Step 4, Temporal Discretization: The Transient Term
------------------------------------------------------

After spatial discretization, the :eq:`eq:fvm_volume_int` can be expressed as,

.. math::
   :label: eq:fvm_diffusionEq_spatial_dis

   \frac{\partial T}{\partial t} V_C - L(T_C^t) =0

where :math:`V_C` is the volume of the discretization cell and :math:`L(T_C^t)` is the spatial discretization operator expressed at some reference time :math:`t`, which can be written as algebraic form (see also :eq:`eq:fvm_matrix_form0`),

.. math::
   :label: eq:fvm_spatial_dis_algebraic

   L(T_C^t) = a_C T_C^t + \sum\limits_{F\sim NB(C)} a_F T_F^t + \sum\limits_{F\sim NB(C)} c_F

where :math:`a_C` is the diagonal coefficients of the matrix, :math:`a_F` is the off-diagonal coefficients, and the :math:`c_F` is the source coefficients as right hand side of matrix of system. 

.. Tip::

   For specific cell :math:`C`, 

   * :math:`a_F` is only contributed from internal faces. 

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


Step 6: Solution of the discretized equations
---------------------------------------------------

The discretization of the differential equation results in a set of discrete algebraic equations, which must be solved to obtain the discrete values of T. The coefficients of these equations may be independent of T (i.e., linear) or dependent on T (i.e. non-linear). The techniques to solve this algebraic system of equations are independent of the discretization method.
The solution methods for solving systems of algebraic equations may be broadly
classified as direct or iterative.

OpenFOAM implementation
===========================

Ok, now let's do some practical tings, 
(0) **generate mesh**; 
(1) **read mesh and do some useful calculation**, e.g. cell volume, face area, ..., etc; 
(2) **discretize Laplacian term and get coefficients**; 
(3) **discretize transient term and get additional coefficients**; 
(4) **construct the final coefficient matrix and RHS**; 
(5) **solve the system of algebraic equations**; 
(6) **write solution to file**.

.. admonition:: Goal

   Deeply look into OpenFOAM and understand how it works!

In order to better understand OpenFOAM's logic and its work flow, we have to look at the basic structure of the `source code <https://cpp.openfoam.org/v6/laplacianFoam_8C_source.html>`_ of a basic solver, :code:`laplacianFoam`.

.. tab:: laplacianFoam.C

   .. code-block::  cpp
      :linenos: 
      :emphasize-lines: 7, 9, 17, 19, 21, 23
      :caption: Source code of laplacianFoam
      :name: lst:source_laplacianFoam

      #include "fvCFD.H"               // Basic head file of OF
      #include "simpleControl.H"       // Basic head file of OF
      int main(int argc, char *argv[]) // Typical c++ main control function
      {
         #include "setRootCaseLists.H" // Do some case/file path-related thing
         #include "createTime.H"       // Create a time object: read controlDict, ...
         #include "createMesh.H"       // (2) Create mesh object: read mesh and do some useful calculation
         simpleControl simple(mesh);   // Time loop control object
         #include "createFields.H"     // (3) Read input data: T and D in the PDE 

         Info<< "\nCalculating temperature distribution\n" << endl;
         while (simple.loop(runTime)) // do time loop
         {
            Info<< "Time = " << runTime.timeName() << nl << endl;
            while (simple.correctNonOrthogonal())  // Non-orthogonal correction loop for the unstructured/non-orthogonal mesh
            {
               fvScalarMatrix TEqn                 // (5) construct the final coefficient matrix, include RHS (.source)
                  (
                     fvm::ddt(T)                   // (4) Discretization of the transient term, return a fvMatrix object
                     - 
                     fvm::laplacian(DT, T)         // (3) Discretization of the Laplacian, return a fvMatrix object
                  );
                  TEqn.solve();                    // (6) Solve the system of algebraic equations
            }
            #include "write.H"                     // (7) Save solution to file.
            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                  << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                  << nl << endl;
         }
         Info<< "End\n" << endl;
         return 0;
      }

   .. list-table:: Comparison between equation and code
      :header-rows: 0

      *  - Term
         - Equation
         - Code
         - Coefficients(ldu:b)
         - Reference
      *  - Transient
         - :math:`\frac{\partial T}{\partial t}`
         - :code:`fvm::ddt(T)`
         - d(:math:`FluxC`), b(:math:`FluxC^o`)
         - :eq:`eq:fvm_firstOrder_Euler_exp_coef`
      *  - Laplacian
         - :math:`\nabla \cdot D \nabla T`
         - :code:`fvm::laplacian(DT, T)`
         - d(:math:`\sum\limits_{F\sim NB(C)}a_F`), u(:math:`a_F` only internal faces)
         - :eq:`eq:fvm_matrix_form_internalCell`, :eq:`eq:fvm_laplacian_coeff_boundary_fixedvalue`, :eq:`eq:fvm_laplacian_coeff_boundary_fixedflux`

.. tab:: createFields.H

   .. code-block::  cpp
      :linenos:
      :emphasize-lines: 0
      :caption: createFields.H
      :name: lst:createFields

      volScalarField T //Read input data of T, include ICs and BCs
      (
         IOobject
         (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
         ),
         mesh
      );

      IOdictionary transportProperties // Read dictionary file of transportProperties
      (
         IOobject
         (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
         )
      );
      dimensionedScalar DT // Read the diffusivity D (constant) from transportProperties object
      (
         transportProperties.lookup("DT")
      );

.. tab:: Input data 1

   .. code-block::  foam
      :linenos:
      :emphasize-lines: 0
      :caption: constant/transportProperties
      :name: lst:transportProperties

      FoamFile
      {
         version     2.0;
         format      ascii;
         class       dictionary;
         location    "constant";
         object      transportProperties;
      }
      DT              DT [0 2 -1 0 0 0 0] 4e-05;

.. tab:: Input data 2

   .. code-block::  foam
      :linenos:
      :emphasize-lines: 0
      :caption: 0/T
      :name: lst:0_T

      FoamFile
      {
         version     2.0;
         format      ascii;
         class       volScalarField;
         object      T;
      }
      dimensions      [0 0 0 1 0 0 0];
      internalField   uniform 273;
      boundaryField
      {
         left
         {
            type            fixedValue;
            value           uniform 273;
         }
         right
         {
            type            fixedValue;
            value           uniform 573;
         }
         "(top|bottom)"
         {
            type            zeroGradient;
         }
      }


Step0, Generate mesh
---------------------------

.. tab:: Generate mesh

   Nothing special, just a meshing process. 
   Here we use the OpenFOAM utility :code:`blockMesh` to generate a regular mesh (:download:`Regular box case <cases/regularBox.zip>`) same as :numref:`fig:polyMesh_regularBox`.

.. tab:: OpenFOAM script

   .. code-block::  foam
      :linenos:
      :emphasize-lines: 0
      :caption: blockMeshDict of a regular box mesh we shown above
      :name: lst:blockMesh_regularBox

      /*--------------------------------*- C++ -*----------------------------------*\
      | =========                 |                                                 |
      | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
      |  \\    /   O peration     | Version:  5                                     |
      |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
      |    \\/     M anipulation  |                                                 |
      \*---------------------------------------------------------------------------*/
      FoamFile
      {
         version     2.0;
         format      ascii;
         class       dictionary;
         object      blockMeshDict;
      }
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

      convertToMeters 1;
      xmin -0.05;    //variable definition
      xmax 0.05;
      ymin -0.015;
      ymax 0.015;
      Lz 0.01;
      vertices    //vertices definition
      (
         ($xmin      $ymin   0)  //coordinate of vertex 0
         ($xmax    $ymin   0)  //coordinate of vertex 1
         ($xmax    $ymax   0)  //coordinate of vertex 2
         ($xmin      $ymax   0)  //coordinate of vertex 3
         ($xmin      $ymin   $Lz)//coordinate of vertex 4
         ($xmax    $ymin   $Lz)//coordinate of vertex 5
         ($xmax    $ymax   $Lz)//coordinate of vertex 6
         ($xmin      $ymax   $Lz)//coordinate of vertex 7
      );
      blocks
      (
         hex (0 1 2 3 4 5 6 7) (10 3 1) simpleGrading (1 1 1)
      );
      boundary
      (
         left    //patch name
         {
            type patch ;
            faces   //face list
            (
                  (0 4 7 3)
            );
         }
         right
         {
            type patch;
            faces
            (
                  (2 6 5 1)
            );
         }
         top
         {
            type patch;
            faces
            (
                  (3 7 6 2)
            );
         }
         bottom
         {
            type patch;
            faces
            (
                  (1 5 4 0)
            );
         }
         frontAndBack    //patch name
         {
            type empty; 
            faces //face list
            (
                  (0 3 2 1)   //back face
                  (4 5 6 7)   //front face
            );
         }
      );
      // ************************************************************************* //

Ok, let's exploring the main steps of the :code:`laplacianFoam`.
Please download (:download:`test_laplacianFoam <cases/test_laplacianFoam.zip>`) and do the following practice/debug steps in :code:`test_laplacianFoam.C`.

Step 1, Read mesh and input field
-------------------------------------

.. tab:: Read data

   The basic properties of cells and faces, e.g. area, volume, face normal vector, will be evaluated after mesh reading, all these processes are happened in the mesh object. 
   It means that after calling :code:`creatFields.H` all these properties are evaluated and stored in the :code:`mesh` object. Of course the temperature field object :code:`T` (volScalarField) with BCs and ICs, and thermal diffusivity :code:`DT` are also initialized from input data after calling :code:`creatFields.H`.
   **It should be noted** that the part of Laplacian discretization coefficients are also calculated in this step. If the mesh is not changed during simulation time, the mesh related coefficients just need to be calculated one time.

   .. list-table:: Mesh- and Field-related coefficients
      :header-rows: 0

      *  - Object
         - Equation
         - Code
         - Faces
         - Reference
      *  - :code:`mesh`
         - :math:`1/\delta_{C\leftrightarrow F}`
         - :code:`mesh.deltaCoeffs()`
         - All faces
         - :eq:`eq:fvm_matrix_form_internalCell`, :eq:`eq:fvm_laplacian_coeff_boundary_fixedvalue`, :eq:`eq:fvm_laplacian_coeff_boundary_fixedflux`
      *  - :code:`T`
         - Dependents on BCs type
         - :code:`gradientInternalCoeffs` (diagonal), :code:`gradientBoundaryCoeffs` (source)
         - Boundary patch/faces
         - :eq:`eq:fvm_laplacian_coeff_boundary_fixedvalue`, :eq:`eq:fvm_laplacian_coeff_boundary_fixedflux`

.. tab:: Access mesh and field properties

   .. code-block::  cpp
      :linenos:
      :emphasize-lines: 0
      :caption: Access mesh properties and field boundary properties
      :name: lst:access_mesh

      // (1). read data
      #include "createMesh.H"
      #include "createFields.H"
      // 1.1 access internal faces
      Info<<"\n\nAccess internal mesh"<<endl;
      surfaceVectorField Cf = mesh.Cf();
      surfaceVectorField Sf = mesh.Sf();
      surfaceScalarField S = mesh.magSf();
      forAll(Cf, iface)
      {
         Info<<iface<<": face center "<<Cf[iface]<<endl;
         Info<<iface<<": face area vector "<<Sf[iface]<<endl;
         Info<<iface<<": face area "<<S[iface]<<endl;
         Info<<iface<<": face delta coeff "<<mesh.deltaCoeffs()[iface]<<endl;
         Info<<iface<<": coeff(D*magSf*deltacoeff) "<<mesh.deltaCoeffs()[iface]*DT*S[iface]<<endl;
      }
      // 1.2. access boundary mesh
      Info<<"\n\nAccess boundary mesh"<<endl;
      const fvBoundaryMesh& boundaryMesh = mesh.boundary(); 
      forAll(boundaryMesh, patchI)
      {
         const fvPatch& patch = boundaryMesh[patchI];
         forAll(patch, faceI)
         {
            Info<<"Patch "<<patch.name()<<" face "<<faceI<<": face center "<<patch.Cf()[faceI]<<endl;
            Info<<"Patch "<<patch.name()<<" face "<<faceI<<": face area vector "<<patch.Sf()[faceI]<<endl;
            Info<<"Patch "<<patch.name()<<" face "<<faceI<<": face area "<<patch.magSf()[faceI]<<endl;
            Info<<"Patch "<<patch.name()<<" face "<<faceI<<": face delta coeff "<<patch.deltaCoeffs()[faceI]<<endl;
            Info<<"Patch "<<patch.name()<<" face "<<faceI<<": owner cell "<<patch.faceCells()[faceI]<<endl;
         } 
      }
      // 1.3. access boundary field, boundary field coefficients, 
      forAll(T.boundaryField(), patchI)
      {
         Info<<"Boundary patch: "<<mesh.boundary()[patchI].name()<<endl;
         Info<<"Is coupled ? "<<mesh.boundary()[patchI].coupled()<<endl;
         Info<<"gradientInternalCoeffs of field T "<<endl;
         Info<<T.boundaryField()[patchI].gradientInternalCoeffs()<<endl; //Diagonal coeff [A]
         Info<<"gradientBoundaryCoeffs of field T "<<endl;
         Info<<T.boundaryField()[patchI].gradientBoundaryCoeffs()<<"\n"<<endl; //source coeff, [B]
      }

.. tab:: Internal mesh

   .. literalinclude::  /_static/log_case/L_FVM/regularMesh_internal.txt
      :language: bash
      :linenos:
      :lines: 0-
      :emphasize-lines: 0
      :caption: Internal mesh properties of the regular mesh shown in :numref:`fig:polyMesh_regularBox`.
      :name: lst:log_internalMesh

.. tab:: Boundary mesh

   .. literalinclude::  /_static/log_case/L_FVM/regularMesh_boundary.txt
      :language: bash
      :linenos:
      :lines: 0-
      :emphasize-lines: 0
      :caption: Boundary mesh properties of the regular mesh shown in :numref:`fig:polyMesh_regularBox`.
      :name: lst:log_boundaryMesh

.. tab:: Boundary T

   .. literalinclude::  /_static/log_case/L_FVM/regularMesh_T_boundary.txt
      :language: bash
      :linenos:
      :lines: 0-
      :emphasize-lines: 0
      :caption: Boundary properties of field T of the regular mesh shown in :numref:`fig:polyMesh_regularBox`.
      :name: lst:log_boundaryT

.. tab:: Internal cell
   :new-set:

   .. figure:: /_figures/Coordinate_delta_internalcell_regularBox.*
      :align: center
      :width: 100 %
      :name: fig:deltaCoeff_InternalCell

      Information of internal cell (:math:`C_{12}`)

.. tab:: Boundary cell

   .. figure:: /_figures/Coordinate_delta_boundary_regularBox.*
      :align: center
      :width: 100 %
      :name: fig:deltaCoeff_BoundaryCell

      Information of boundary cell (:math:`C_{19}`)

.. _OF_fvmLaplacian:

Step 2, discretize Laplacian term
-------------------------------------

.. tab:: discretize Laplacian term

   Because discretization coefficients matrix of Laplacian term is a symmetric matrix, so :code:`fvm::Laplacian(DT, T)` will return a fvMatrix object only has diagonal and upper. 
   What :code:`fvm::Laplacian` actually did is evaluate (1) :math:`a_F` (see :eq:`eq:fvm_matrix_form_internalCell`) for each internal faces, (2) :math:`a_C` for each cells, which is the negative summation of :math:`a_F`, (3) store the field BCs-related coefficients as :code:`internalCoeffs` and :code:`boundaryCoeffs`, respectively.
   All of these are implemented in gaussLaplacianScheme.C_ . **Note that** the :code:`Gaussian` scheme is the only choice for Laplacian discretization in OF. 

.. tab:: Access Laplacian coeffs 

   .. code-block::  cpp
      :linenos:
      :emphasize-lines: 0
      :caption: Access fvm::Laplacian discretization
      :name: lst:access_Laplacian

      fvScalarMatrix Laplacian(fvm::laplacian(DT, T));
      Info<<"fvm::laplacian(DT, T): "<<"\n"
         <<"\tLower"<<Laplacian.lower()<<"\n"
         <<"\tDiagonal"<<Laplacian.diag()<<"\n"
         <<"\tUpper"<<Laplacian.upper()<<"\n"
         <<"\tinternalCoeffs"<<Laplacian.internalCoeffs()<<"\n"
         <<"\tboundaryCoeffs"<<Laplacian.boundaryCoeffs()<<"\n"
         <<"\tSource"<<Laplacian.source()<<"\n"
         <<endl;

   .. tip::

      The source coefficients come from BCs are not stored in the :code:`.source()`, but in :code:`boundaryCoeffs`.
      So if you print :code:`Laplacian.source()`, it will display zero.
      The BCs-related source will be added into the :code:`TEqn.source` when :code:`TEqn.solve()` is calling. 
      There is a protected member function named :code:`addBoundarySource` will be called in :code:`solve()` function.

.. tab:: Source code of Laplacian

   .. code-block::  cpp
      :linenos:
      :emphasize-lines: 15, 20, 39, 40
      :caption: Key source code of fvm::Laplacian (gaussLaplacianScheme.C_ before nonOrthogonal correcting)
      :name: lst:source_fvm_Laplacian

      template<class Type, class GType>
      tmp<fvMatrix<Type>>
      gaussLaplacianScheme<Type, GType>::fvmLaplacianUncorrected
      (
         const surfaceScalarField& gammaMagSf,
         const surfaceScalarField& deltaCoeffs,
         const GeometricField<Type, fvPatchField, volMesh>& vf
      )
      {
         tmp<fvMatrix<Type>> tfvm
         (
            new fvMatrix<Type>
            (
                  vf,
                  deltaCoeffs.dimensions()*gammaMagSf.dimensions()*vf.dimensions()
            )
         );
         fvMatrix<Type>& fvm = tfvm.ref();

         fvm.upper() = deltaCoeffs.primitiveField()*gammaMagSf.primitiveField();
         fvm.negSumDiag();

         forAll(vf.boundaryField(), patchi)
         {
            const fvPatchField<Type>& pvf = vf.boundaryField()[patchi];
            const fvsPatchScalarField& pGamma = gammaMagSf.boundaryField()[patchi];
            const fvsPatchScalarField& pDeltaCoeffs =
                  deltaCoeffs.boundaryField()[patchi];

            if (pvf.coupled())
            {
                  fvm.internalCoeffs()[patchi] =
                     pGamma*pvf.gradientInternalCoeffs(pDeltaCoeffs);
                  fvm.boundaryCoeffs()[patchi] =
                     -pGamma*pvf.gradientBoundaryCoeffs(pDeltaCoeffs);
            }
            else
            {
                  fvm.internalCoeffs()[patchi] = pGamma*pvf.gradientInternalCoeffs();
                  fvm.boundaryCoeffs()[patchi] = -pGamma*pvf.gradientBoundaryCoeffs();
            }
         }

         return tfvm;
      }

.. tab:: Laplacian coefficients

   .. literalinclude::  /_static/log_case/L_FVM/regularMesh_fvm_Laplacian.txt
      :language: bash
      :linenos:
      :lines: 0-
      :emphasize-lines: 58, 63, 67, 68, 77, 70
      :caption: fvm::Laplacian Coefficients  of the regular mesh shown in :numref:`fig:polyMesh_regularBox`.
      :name: lst:log_fvm_Laplacian

.. tab:: Internal cell
   :new-set:

   .. figure:: /_figures/Coordinate_Laplacian_internalcell_regularBox.*
      :align: center
      :width: 100 %
      :name: fig:LaplacianCoeff_InternalCell

      Information of internal cell (:math:`C_{12}`)

.. tab:: Boundary cell 19

   .. figure:: /_figures/Coordinate_Laplacian_boundary_C19_regularBox.*
      :align: center
      :width: 100 %
      :name: fig:LaplacianCoeff_BoundaryCell

      Information of boundary cell (:math:`C_{19}`)

.. tab:: Boundary cell 9

   .. figure:: /_figures/Coordinate_Laplacian_boundary_C9_regularBox.*
      :align: center
      :width: 100 %

      Information of boundary cell (:math:`C_{9}`)

.. tab:: Boundary cell 0

   .. figure:: /_figures/Coordinate_Laplacian_boundary_C0_regularBox.*
      :align: center
      :width: 100 %

      Information of boundary cell (:math:`C_{0}`)

.. tab:: Boundary cell 10

   .. figure:: /_figures/Coordinate_Laplacian_boundary_C10_regularBox.*
      :align: center
      :width: 100 %

      Information of boundary cell (:math:`C_{10}`)

.. tab:: Boundary cell 5

   .. figure:: /_figures/Coordinate_Laplacian_boundary_C5_regularBox.*
      :align: center
      :width: 100 %

      Information of boundary cell (:math:`C_{5}`)

.. _OF_fvmDdt:

Step 3, discretize transient term
-------------------------------------



.. tab:: discretize Laplacian term

   For implicit discretization, :code:`fvm::ddt(T)` will return a fvMatrix object contains diagonal coefficients and source.
   The coefficients depend on discretization scheme.
   For example Euler scheme, the diagonal coefficients are calculated from :eq:`eq:fvm_temporal_linear` and :eq:`eq:fvm_firstOrder_Euler_imp_coeff`.
   All these are implemented in EulerDdtScheme.C_. 

   .. tip::

      There are 8 schemes for transient discretization in OpenFOAM

      #. CoEuler
      #. CrankNicolson
      #. Euler
      #. SLTS
      #. backward
      #. bounded
      #. localEuler
      #. steadyState

.. tab:: Access fvm::ddt coeffs

   .. code-block::  cpp
      :linenos:
      :emphasize-lines: 0
      :caption: Access fvm::ddt discretization
      :name: lst:access_ddt

      Info<<"fvm::ddt(T): "<<"\n"
         <<"\tLower"<<ddt.lower()<<"\n"
         <<"\tDiagonal"<<ddt.diag()<<"\n"
         <<"\tUpper"<<ddt.upper()<<"\n"
         <<"\tinternalCoeffs"<<ddt.internalCoeffs()<<"\n" //actually this is not necessary for fvm::ddt, this is definitely equal to zero
         <<"\tboundaryCoeffs"<<ddt.boundaryCoeffs()<<"\n"
         <<"\tSource"<<ddt.source()<<"\n"
         <<endl;

.. tab:: Source code of fvm::ddt

   .. code-block::  cpp
      :linenos:
      :emphasize-lines: 19, 21, 29
      :caption: Key source code of fvm::ddt (EulerDdtScheme.C_ )
      :name: lst:source_fvm_ddt

      template<class Type>
      tmp<fvMatrix<Type>>
      EulerDdtScheme<Type>::fvmDdt
      (
         const GeometricField<Type, fvPatchField, volMesh>& vf
      )
      {
         tmp<fvMatrix<Type>> tfvm
         (
            new fvMatrix<Type>
            (
                  vf,
                  vf.dimensions()*dimVol/dimTime
            )
         );

         fvMatrix<Type>& fvm = tfvm.ref();

         scalar rDeltaT = 1.0/mesh().time().deltaTValue(); // 1/dt

         fvm.diag() = rDeltaT*mesh().Vsc(); // Vc/dt (FluxC)

         if (mesh().moving())
         {
            fvm.source() = rDeltaT*vf.oldTime().primitiveField()*mesh().Vsc0();
         }
         else
         {
            fvm.source() = rDeltaT*vf.oldTime().primitiveField()*mesh().Vsc(); // -T_old*Vc/dt
         }

         return tfvm;
      }

.. tab:: fvm::ddt coefficients

   :math:`\Delta t = 0.05\ s`

   .. literalinclude::  /_static/log_case/L_FVM/regularMesh_fvm_ddt.txt
      :language: bash
      :linenos:
      :lines: 0-
      :emphasize-lines: 0
      :caption: fvm::ddt Coefficients  of the regular mesh shown in :numref:`fig:polyMesh_regularBox`.
      :name: lst:log_fvm_ddt


Step 4, construct the final coefficient matrix and RHS
-----------------------------------------------------------

.. tab:: Final matrix

   The final coefficient matrix is constructed by simply adding the matrix of :ref:`OF_fvmLaplacian` and :ref:`OF_fvmDdt`.
   
   * The diagonal coefficients come from :math:`a_C` of Laplacian term and :math:`FluxC` of transient term
   * The off-diagonal coefficients only come from :math:`a_F (internal\ face)` of Laplacian term
   * The RHS comes from :math:`c_F` of Laplacian term when at boundary faces and :math:`FluxC^oT_C^o` of transient term.

.. tab:: Coefficients at the first time step

   .. literalinclude::  /_static/log_case/L_FVM/regularMesh_fvm_TEqn.txt
      :language: bash
      :linenos:
      :lines: 0-
      :emphasize-lines: 69, 76, 57, 66, 62, 67
      :caption: Final matrix coefficients of the regular mesh shown in :numref:`fig:polyMesh_regularBox` at the first time step.
      :name: lst:log_fvm_TEqn

Step 5, solve
-------------------------------------

For the :download:`Regular box case <cases/regularBox.zip>` case, we can use **PBiCG** solver and **DILU** preconditioner.

.. admonition:: Available preconditioner in OpenFOAM

   * **diagonal** : for symmetric & nonsymmetric matrices (not very effective)
   * **DIC** : Diagonal Incomplete Cholesky preconditioner for symmetric matrices
   * **DILU** : Diagonal Incomplete LU preconditioner for nonsymmetric matrices
   * **FDIC** : Fast Diagonal Incomplete Cholesky preconditioner
   * **GAMG** : Geometric Agglomerated algebraic MultiGrid preconditioner


.. admonition:: Available solver in OpenFOAM

   * **BICCG**: Diagonal incomplete LU preconditioned BiCG solver
   * **diagonalSolver**: Solver for symmetric and nonsymmetric matrices
   * **GAMG**: Geometric Agglomerated algebraic Multi-Grid solver
   * **ICC**: Incomplete Cholesky Conjugate Gradient solver
   * **PBiCG**: Bi-Conjugate Gradient solver with preconditioner
   * **PCG**: Conjugate Gradient solver with preconditioner
   * **smoothSolver**: Iterative solver with run-time selectable smoother

.. admonition:: Krylov-subspace solvers

   * **CG**: The Conjugate Gradient algorithm applies to systems where A is symmetric positive definite (SPD)
   * **GMRES**: The Generalized Minimal RESidual algorithm is the first method to try if A is not SPD.
   * **BiCG**: The BiConjugate Gradient algorithm applies to general linear systems, but the convergence can be quite erratic.
   * **BiCGstab**: The stabilized version of the BiConjugate Gradient algorithm.



Step 6, write
-------------------------------------




Jupyter notebook
-------------------

.. toctree::
    :maxdepth: 2

    cases/jupyter/VisualizeResults.ipynb


