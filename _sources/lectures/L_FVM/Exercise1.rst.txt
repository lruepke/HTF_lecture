.. include:: /include.rst_

.. _L_FVM_Exercise1:

Exercise 1
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
    
    b_C = \iint\limits_{S_C} (D\nabla T)\cdot d\vec{S}

The :eq:`eq:fvm_surface_int` is actually a heat balance over cell :math:`C`. 
It is basically the integral form of the original partial differential equation and involves **no approximation**.

3. Transform the surface integral in :eq:`eq:fvm_surface_int` as a summation over the control volume faces (**here introduce the first approximation**),

.. math::
    :label: eq:fvm_surface_sum
    
    b_C = \sum\limits_{f\sim nb(C)} (D\nabla T)_f\cdot \vec{S}_f

Here :math:`f` denotes the integration point at the centroid of the bounding face.
For a specific cell :math:`C` shown in :numref:`fvm_geometry_bcs_regularBox`, 
the :eq:`eq:fvm_surface_sum` can be expanded as,

.. math::
    :label: eq:fvm_surface_sum_expand
    
    b_C = \color{red}{(D\nabla T)_{f_1} \cdot \vec{S}_{f_1}} + (D\nabla T)_{f_2} \cdot \vec{S}_{f_2} + (D\nabla T)_{f_3} \cdot \vec{S}_{f_3}+ (D\nabla T)_{f_4} \cdot \vec{S}_{f_4}

4. Evaluating gradient at centroid of faces

   4.1 Considering face :math:`f_1`, calculate the first term on the right hand side of the :eq:`eq:fvm_surface_sum_expand`,

   .. math::
      :label: eq:fvm_surface_term_f1
      
      \color{red}{(D\nabla T)_{f_1} \cdot \vec{S}_{f_1}} = \left(D \frac{\partial T}{\partial x} \vec{i} + D\frac{\partial T}{\partial y} \vec{j} \right)_{f_1} \cdot \Delta y_{f_1} \vec{i} = \color{blue}{\left( D\frac{\partial T}{\partial x} \right)_{f_1} \Delta y_{f_1}}

   where 

   * :math:`S_{f_1} = \Delta y_{f_1}`  is the area of face  :math:`f_1` by assuming :math:`\Delta z = 1`. 
   * :math:`\vec{S}_{f_1} = S_{f_1} \vec{n}_{f_1}` is the surface vector of face :math:`f_1`.
   * :math:`\vec{n}_{f_1} = \vec{i}` is the normal vector of the face :math:`f_1` directed out of the cell :math:`C` (see :numref:`fig:polyMesh_regularBox`).

   4.2 Calculate gradient of :math:`T` at face centroid, here introduce **the second approximation**, e.g. assuming linear variation of T and then the gradient term in the :eq:`eq:fvm_surface_term_f1` can be approximated as,

   .. math::

      \color{blue}{\left( D\frac{\partial T}{\partial x} \right)_{f_1} \Delta y_{f_1}} = D\frac{T_{F_1} - T_{C}}{x_{F_1} - x_C}\Delta y_{f_1} = D\frac{T_{F_1} - T_{C}}{\delta x_{f_1}}\Delta y_{f_1} = \color{green}{a_{F_1}} (T_{F_1} - T_C)

   where :math:`\delta x_{f_1}` represents the distance between cell :math:`C` and :math:`F_1`.

   4.3 Do the same thing for the remain faces (:math:`f_2, f_3, f_4`) and get their coefficients :math:`\color{green}{a_{F_2}, a_{F_3}, a_{F_4}}` respectively,

   .. math::
      :label: eq:fvm_aFaC

      \color{green}{a_{F_1}} = \frac{\Delta y_{f_1}}{\delta x_{f_1}},
      \color{green}{a_{F_2}} = \frac{\Delta x_{f_2}}{\delta y_{f_2}},
      \color{green}{a_{F_3}} = \frac{\Delta x_{f_3}}{\delta y_{f_3}},
      \color{green}{a_{F_4}} = \frac{\Delta y_{f_4}}{\delta x_{f_4}}

5. Construct matrix form of the :eq:`eq:fvm_surface_sum` by substituting the above relations,

.. math::
   :label: eq:fvm_matrix_form0

   \begin{eqnarray}
   \sum\limits_{f\sim nb(C)} (D\nabla T)_f\cdot \vec{S}_f & = & \sum\limits_{F\sim NB(C)} a_F(T_F - T_C) = \sum\limits_{F\sim NB(C)} a_F T_F - T_C\sum\limits_{F\sim NB(C)} a_F \\
   & = & a_CT_C + \sum\limits_{F\sim NB(C)} a_F T_F
   \end{eqnarray}

Considering :math:`b_C` matrix form of the :eq:`eq:fvm_surface_sum` can be written as,

.. math::
   :label: fvm_matrix_form

   \mathbf{A}_{N \times N} \mathbf{T}_{N\times1} = \mathbf{b}_{N\times1}

here :math:`N` is the number of cell, :math:`\mathbf{A}` is a sparse matrix of coefficients, 
:math:`\mathbf{T}` is cell-centered temperature field and :math:`\mathbf{b}` is the transient term which will explained in the later section.
For the steady state problem, :math:`\mathbf{b} \equiv 0`. 
The matrix is visualized as :numref:`fig:fvm_matrix`.

.. tab:: Regular mesh 

   .. figure:: /_figures/matrix_FVM_regularBox.*
      :align: center
      :name: fig:fvm_matrix
      :figwidth: 100%

      Visualization of :math:`\mathbf{A}_{N \times N} \mathbf{T}_{N\times1} = \mathbf{b}_{N\times1}`. 
      The coefficients of the selected cell :math:`C` (:numref:`fig:polyMesh_regularBox`) are marked by green rectangle.

.. tab:: unstructured mesh 

   .. figure:: /_figures/matrix_FVM_unstructured.*
      :align: center
      :name: fig:fvm_matrix_unstructured
      :figwidth: 100%

      Visualization of :math:`\mathbf{A}_{N \times N} \mathbf{T}_{N\times1} = \mathbf{b}_{N\times1}`. 
      The coefficients of the selected cell :math:`C` (:numref:`fig:polyMesh_regularBox`) are marked by green rectangle.

Jupyter notebook
-------------------

.. toctree::
    :maxdepth: 2

    cases/jupyter/VisualizeResults.ipynb


