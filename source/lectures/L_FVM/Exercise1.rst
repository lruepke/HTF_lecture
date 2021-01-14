.. include:: /include.rst_

.. _L_FVM_Exercise1:

.. warning::

   Before we start the exercise 1, please download (:download:`test_laplacianFoam <cases/test_laplacianFoam.zip>`) the modified Laplacian solver first.
   Then put it in the shared folder and compile (:code:`wmake`) it in the docker container.

==============
Exercise 1
==============


Theory
===========


.. admonition:: Goals

   - explore FVM theory
   - link FVM theory and OpenFOAM solver implementation
   - step by step exploring details of the solver and solving process


Step 1: Geometric and physical modeling
---------------------------------------------------

.. math::
    :label: eq:fvm_laplacian_dif
    
    \frac{\partial T}{\partial t} = \nabla \cdot D \nabla T 

Here :math:`D` is the thermal diffusivity, can be a constant value.
see also :ref:`theory_heat_diffusion`

.. figure:: /_figures/boundaryConditions_FVM_regularBox.*
   :align: center
   :name: fvm_geometry_bcs_regularBox
   :figwidth: 100%

   Model geometry, boundary conditions and initial condition.

.. admonition:: What we need to solve ?

   Temperature field :math:`T` of the model region at time :math:`t`, 
   therefore :math:`T` is the only unknown variable of our problem.


Step 2: Domain discretization
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
   The normal vector of face always points from the **owner** cell to the **neighbour** cell.

Step 3: Equation discretization
---------------------------------------------------

The equation discretization step is performed over each element of the computational domain to yield an algebraic relation that connects the value of a variable in an element to the values of the variable in the neighboring elements :cite:`moukalled2016finite`. 

.. admonition:: Goal

   Transform partial differential :eq:`eq:fvm_laplacian_dif` into a set of algebraic equations: :math:`\mathbf{A}[T] = \mathbf{b}`.

1. Integrate :eq:`eq:fvm_laplacian_dif` over element :math:`C` (the green cell in :numref:`fig:polyMesh_regularBox`) that enables recovering its integral balance form as,

.. math::
    :label: eq:fvm_volume_int
    
    \iiint\limits_{V_C} \frac{\partial T}{\partial t} dV \equiv b_C  = \iiint\limits_{V_C} \nabla \cdot D\nabla T dV

**Note** that the transient term :math:`\iiint\limits_{V_C} \frac{\partial T}{\partial t} dV \equiv b_C` will be processed in the later section

2. Transform the volume integral on the right hand side of :eq:`eq:fvm_volume_int` into a surface integral by applying the divergence theorem, 

.. math::
    :label: eq:fvm_surface_int
    
    b_C = \oint\limits_{\partial V_C} (D\nabla T)\cdot d\vec{S}

The :eq:`eq:fvm_surface_int` is actually a heat balance over cell :math:`C`. 
It is basically the integral form of the original partial differential equation and involves **no approximation**.

Here we can introduce definition of **heat diffusion flux**, 

.. math::
   :label: eq:fvm_flux_D

   \vec{J}^{T,D} \equiv D\nabla T 

3. Transform the surface integral in :eq:`eq:fvm_surface_int` as a summation over the control volume faces (**still no approximation**),

.. math::
    :label: eq:fvm_surface_sum
    
    b_C = \sum\limits_{f\sim faces(V_C)} \iint\limits_{f}\vec{J}^{T,D}_f \cdot \vec{S}_f 

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

      b_{C_{12}} = \color{red}{\iint_{f_{23}}\vec{J}^{T,D}_{f_{23}} \cdot \vec{S}_{f_{23}}} + \iint_{f_{24}}\vec{J}^{T,D}_{f_{24}} \cdot \vec{S}_{f_{24}} + \iint_{f_{21}}\vec{J}^{T,D}_{f_{21}} \cdot \vec{S}_{f_{21}}  + \iint_{f_{5}}\vec{J}^{T,D}_{f_{5}} \cdot \vec{S}_{f_{5}} 


   4.1 Considering face :math:`f_{23}`, calculate the first term on the right hand side of the :eq:`eq:fvm_surface_sum_expand` (we have to introduce the **first approximation** at this step). Using a Gaussian quadrature the integral at the face :math:`f_{23}`, for example, becomes,

   .. math::
      :label: fvm_flux_Gaussian_integral

      \color{red}{\iint_{f_{23}}\vec{J}^{T,D}_{f_{23}} \cdot \vec{S}_{f_{23}}} = \iint_{f_{23}} (\vec{J}^{T,D}_{f_{23}} \cdot \vec{n}_{f_{23}}) dS_{f_{23}} \approx \color{orange}{\sum\limits_{ip\sim ip(f_{23})} (\vec{J}^{T,D} \cdot \vec{n}_{ip})\omega_{ip} S_{f_{23}}}

   where :math:`S_f` is the area of face :math:`f_{23}`, :math:`ip` refers to a integration point and :math:`ip(f_{23})` the number of integration points along surface :math:`f_{23}`, :math:`\omega_{ip}` is the integral weights. 

   4.2 Choose integral scheme or integral points

   To simply explain the calculation process and logic, here we adopt an **one integration point** scheme with weight :math:`\omega = 1`, thus :eq:`fvm_flux_Gaussian_integral` becomes,

   .. math::
      :label: eq:fvm_surface_term_f1
      
      \color{orange}{\sum\limits_{ip\sim ip(f_{23})} (\vec{J}^{T,D} \cdot \vec{n}_{ip})\omega S_f } = \left(D \frac{\partial T}{\partial x} \vec{i} + D\frac{\partial T}{\partial y} \vec{j} \right)_{f_{23}} \cdot \Delta y_{f_{23}} \vec{i} = \color{blue}{\left( D\frac{\partial T}{\partial x} \right)_{f_{23}} \Delta y_{f_{23}}}

   where 

   * :math:`S_{f_{23}} = \Delta y_{f_{23}}`  is the area of face  :math:`f_{23}` by assuming :math:`\Delta z = 1`. 
   * :math:`\vec{S}_{f_{23}} = S_{f_{23}} \vec{n}_{f_{23}}` is the surface vector of face :math:`f_{23}`.
   * :math:`\vec{n}_{f_{23}} = \vec{i}` is the normal vector of the face :math:`f_{23}` directed out of the cell :math:`C` (see :numref:`fig:polyMesh_regularBox`).

   4.3 Calculate gradient of :math:`T` at face centroid, here introduce **the second approximation**, e.g. assuming linear variation of T and then the gradient term in the :eq:`eq:fvm_surface_term_f1` can be approximated as,

   .. math::

      \color{blue}{\left( D\frac{\partial T}{\partial x} \right)_{f_{23}} \Delta y_{f_{23}}} = D\frac{T_{f_{23}} - T_{C}}{x_{f_{23}} - x_C}\Delta y_{f_{23}} = D\frac{T_{f_{23}} - T_{C}}{\delta x_{f_{23}}}\Delta y_{f_{23}} = \color{green}{a_{f_{23}}} (T_{f_{23}} - T_C)

   where :math:`\delta x_{f_{23}}` represents the distance between cell :math:`C` and :math:`f_{23}`.

   Do the same thing for the remain faces (:math:`f_2, f_3, f_4`) and get their coefficients :math:`\color{green}{a_{F_2}, a_{F_3}, a_{F_4}}` respectively,

   .. math::
      :label: eq:fvm_aFaC

      \color{green}{a_{f_{23}}} = D\frac{\Delta y_{f_{23}}}{\delta x_{f_{23}}},
      \color{green}{a_{F_2}} = D\frac{\Delta x_{f_2}}{\delta y_{f_2}},
      \color{green}{a_{F_3}} = D\frac{\Delta x_{f_3}}{\delta y_{f_3}},
      \color{green}{a_{F_4}} = D\frac{\Delta y_{f_4}}{\delta x_{f_4}}
   
   4.3 Construct coefficients of a specific cell :math:`C_{12}` 

   Now we get all the coefficients, thus the :eq:`eq:fvm_surface_sum_expand` can be discretized and expressed as a matrix form,

   .. math::
      :label: eq:fvm_matrix_form_internalCell
      
      \begin{eqnarray}
      b_{C_{12}} & =  & a_{F_{23}} (T_{F_{23}} - T_C) + a_{F_{24}} (T_{F_{24}} - T_C) + a_{F_{21}} (T_{F_{21}} - T_C) + a_{F_{5}} (T_{F_{5}} - T_C)\\
         & = & \sum\limits_{f\in[23, 24, 21, 5]} a_{F_f} T_{F_f} + \left(-\sum\limits_{f\in[23, 24, 21, 5]} a_{F_f} \right) T_C \\
         & = & \sum\limits_{f\sim nb(\text{internal face})} a_{F_f} T_{F_f} + \left(-\sum\limits_{f\sim nb(\text{all face})} a_{F_f} \right) T_C  \\
      \end{eqnarray}
   
   .. math:: 
      :label: eq:fvm_matrix_form_internalCell12

      \mathbf{A}_{N\times N}[12,:] \mathbf{T}_{N\times 1} =  B_{12}

   where :math:`B_{12} = b_{12}`.
   
.. tab:: Boundary cell

   For a specific boundary cell, e.g. **cell 19** shown in :numref:`fig:polyMesh_regularBox`, the :eq:`eq:fvm_surface_sum` can be rewritten as,

   .. math::
      :label: eq:fvm_surface_sum_expand_boundary

      b_{C_{19}} = \iint_{f_{18}}\vec{J}^{T,D}_{f_{18}} \cdot \vec{S}_{f_{18}} + \iint_{f_{35}}\vec{J}^{T,D}_{f_{35}} \cdot \vec{S}_{f_{35}} +\iint_{f_{37}}\vec{J}^{T,D}_{f_{37}} \cdot \vec{S}_{f_{37}} + \color{red}{\iint_{f_{51}}\vec{J}^{T,D}_{f_{51}} \cdot \vec{S}_{f_{51}}}
   
   The first and second terms on the right hand side are flux of internal face, so there is nothing special and just calculate them following steps for internal faces we did before. The special thing is the boundary face :math:`f_{17}`, we can call it :math:`f_b` for a general case. **Now the problem is how to evaluate flux on the boundary face** :math:`f_b`.

   The The fluxes on the interior faces are discretized as before, while **the boundary flux is discretized with the aim of constructing a linearization with respect to** the cell field :math:`T_C`, e.g. :math:`T_{C_{19}}` of cell :math:`C_{19}` shown in :numref:`fig:polyMesh_regularBox`, thus 

   .. math::
      :label: eq:fvm_surface_flux_boundary

      \color{red}{\iint_{f_{b}}\vec{J}^{T,D}_{f_{b}} \cdot \vec{S}_{f_{b}}} = -a_{F_b} T_C + c_{F_b}

   All right, **now out goal is to determine the coefficients of** :math:`a_{F_b}` and :math:`c_{F_b}` according to the boundary conditions. Generally there are two basic kinds of boundary conditions, they are **Dirichlet boundary condition** (also called the first type or fixed value boundary condition) and **Von Neumann boundary condition** (all called the second type or fixed flux boundary condition), respectively.

   4.1 Fixed value boundary condition

   Fixed value means the value of field :math:`T` is given on the boundary face is give, i.e. :math:`\color{purple}{T_{f_b}}` is given. There fore we can calculate flux on boundary face, :math:`f_{51}` for example,

   .. math::
      :label: eq:fvm_flux_boundary

      \color{red}{\iint_{f_{51}}\vec{J}^{T,D}_{f_{51}} \cdot \vec{S}_{f_{51}}} = \iint_{f_{51}} (\vec{J}^{T,D}_{f_{51}} \cdot \vec{n}_{f_{51}}) dS_{f_{51}} \approx \color{orange}{\sum\limits_{ip\sim ip(f_{51})} (\vec{J}^{T,D} \cdot \vec{n}_{ip})\omega_{ip} S_{f_{51}}}

   Here we still use one integration point scheme to explain the calculation process,

   .. math::
      :label: eq:fvm_flux_Gaussian_integral_boundary

      \color{orange}{\sum\limits_{ip\sim ip(f_{51})} (\vec{J}^{T,D} \cdot \vec{n}_{ip})\omega S_f } = \left(D \frac{\partial T}{\partial x} \vec{i} + D\frac{\partial T}{\partial y} \vec{j} \right)_{f_{51}} \cdot \Delta y_{f_{51}} \vec{i} = \color{blue}{\left( D\frac{\partial T}{\partial x} \right)_{f_{51}} \Delta y_{f_{51}}}

   Now let's calculate the gradient on the boundary face,

   .. math::
      :label: eq:fvm_gradient_boundary

      \color{blue}{\left( D\frac{\partial T}{\partial x} \right)_{f_{51}} \Delta y_{f_{51}}} = D\frac{\color{purple}{T_{f_{51}}} - T_{C}}{x_{f_{51}} - x_C}\Delta y_{f_{51}} = D\frac{\color{purple}{T_{f_{51}}} - T_{C}}{\delta x_{f_{51}}}\Delta y_{f_{51}} = \color{green}{a_{F_{51}}} (T_{f_{51}} - T_C)

   Substituting :eq:`eq:fvm_gradient_boundary` back to  :eq:`eq:fvm_flux_Gaussian_integral_boundary`, :eq:`eq:fvm_flux_boundary` and  :eq:`eq:fvm_surface_flux_boundary`, we can get Laplacian discretization coefficients :math:`a_{f_{51}}` and :math:`c_{f_{51}}`,

   .. math::
      :label: eq:fvm_laplacian_coeff_boundary_fixedvalue

      \begin{eqnarray}
      a_{F_{51}} &=& D\frac{\Delta y_{f_{51}}}{\delta x_{f_{51}}}\\
      c_{F_{51}} &=& a_{F_{51}}\color{purple}{T_{f_{51}}}
      \end{eqnarray}

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
      b_{C_{19}} & =  & a_{F_{18}} (T_{F_{18}} - T_C) + a_{F_{35}} (T_{F_{35}} - T_C) + a_{F_{37}} (T_{F_{37}} - T_C) + a_{F_{51}} T_C +c_{F_{51}} \\
         & = & \sum\limits_{f\in[18, 35, 37]} a_{F_f} T_{F_f} + \left(-\sum\limits_{f\in[18, 35, 37, 51]} a_{F_f} \right) T_C + \sum\limits_{f\in[51]} c_{F_f} \\
         & = & \sum\limits_{f\sim nb(\text{internal face})} a_{F_f} T_{F_f} + \left(-\sum\limits_{f\sim nb(\text{all face})} a_{F_f} \right) T_C + \sum\limits_{f\sim nb(\text{boundary face})} c_{F_f}
      \end{eqnarray}

   .. math:: 
      :label: eq:fvm_matrix_form_boundaryCell19

      \mathbf{A}_{N\times N}[19,:] \mathbf{T}_{N\times 1} =  B_{19}

   where :math:`B_{12} = b_{C_{12}} - \sum\limits_{f\sim nb(\text{boundary face})} c_{F_f}`.
      
5. Construct general matrix form of the :eq:`eq:fvm_surface_sum`

.. math::
   :label: eq:fvm_matrix_form0

   \begin{eqnarray}
   \sum\limits_{f\sim nb(C)} (D\nabla T)_{F_f}\cdot \vec{S}_{F_f} & = &  \sum\limits_{f\sim \text{internal faces}(C)} a_{F_f} T_{F_f} + \left(-\sum\limits_{f\sim nb(C)} a_{F_f} \right)T_C + \sum\limits_{f\sim \text{boundary faces}(C)} c_{F_{F_f}} \\
   & = & a_CT_C + \sum\limits_{f\sim \text{internal faces}(C)} a_{F_f} T_{F_f} + \sum\limits_{f\sim \text{boundary faces}(C)} c_{F_{f}}
   \end{eqnarray}

Matrix form of the :eq:`eq:fvm_surface_sum` can be written as,

.. math::
   :label: fvm_matrix_form

   \mathbf{A}_{N \times N} \mathbf{T}_{N\times1} = \mathbf{B}_{N\times1}

here :math:`N` is the number of cell, :math:`\mathbf{A}` is a sparse matrix of coefficients, 
:math:`\mathbf{T}` is cell-centered temperature field.
The matrix is visualized as :numref:`fig:fvm_matrix`.
It should be noted that **the coefficients matrix of Laplacian term is a symmetric matrix and the boundary conditions only affect diagonal entry of the coefficients matrix**. Further, only the fixed value boundary condition affects the matrix and the fixed flux boundary condition affects the RHS :math:`B`.

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


Step 5: Transient discretization 
---------------------------------------------------

.. warning::

   working on this


Step 6: Solution of the discretized equations
---------------------------------------------------

The discretization of the differential equation results in a set of discrete algebraic equations, which must be solved to obtain the discrete values of T. The coefficients of these equations may be independent of T (i.e., linear) or dependent on T (i.e. non-linear). The techniques to solve this algebraic system of equations are independent of the discretization method.
The solution methods for solving systems of algebraic equations may be broadly
classified as direct or iterative.

OpenFOAM implementation
===========================

Step 1: Read data
--------------------------

Read mesh 
^^^^^^^^^^^^^

.. admonition:: 10 seconds thinking

   What information do we need in the mesh object ?

.. #. Surface centroids

.. .. math::
..    :label: fvm_mesh_C

..    (\mathbf{X}_{Centorid})_f = \frac{\sum\limits_{t \sim Sub-triangles(f)} (\mathbf{X}_{centroid})_t S_t}{S_f}
   
.. tab:: FVM mesh

   #. Centroid of face
   #. Area of face
   #. Vector of face
   #. Centroid of cell
   #. Volume of cell

.. tab:: OpenFOAM code snippit

   .. code-block:: cpp

      surfaceVectorField Cf = mesh.Cf();
      forAll(Cf, iface)
      {
         Info<<Cf[iface][0]<<" "<<Cf[iface][1]<<" "<<Cf[iface][2]<<endl;
      }

Jupyter notebook
-------------------

.. toctree::
    :maxdepth: 2

    cases/jupyter/VisualizeResults.ipynb


