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
    
    \frac{\partial T}{\partial t} = D \nabla^2 T 

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

.. figure:: /_figures/mesh_FVM_regularBox.*
   :align: center
   :name: fig:polyMesh_regularBox
   :figwidth: 100%

   2D regular mesh structure and topology information of OpenFOAM polyMesh (See :ref:`/lectures/L_FVM/cases/regularBox/jupyter/VisualizeResults.ipynb#1.-Read-and-plot-mesh-information` in the notebook).

Step 3: Equation discretization
---------------------------------------------------

The equation discretization step is performed over each element of the computational domain to yield an algebraic relation that connects the value of a variable in an element to the values of the variable in the neighboring elements :cite:`moukalled2016finite`. 

.. admonition:: Goal

   Transform partial differential equation :eq:`eq:fvm_laplacian_dif` into a set of algebraic equations: :math:`\mathbf{A}[T] = \mathbf{b}`.

1. Integrate equation :eq:`eq:fvm_laplacian_dif` over element :math:`C` (:numref:`fig:polyMesh_regularBox`) that enables recovering its integral balance form as,

.. math::
    :label: eq:fvm_volume_int
    
    \iint\limits_{V_C} dV = \iint\limits_{V_C} D\nabla^2 T dV

Jupyter notebook
-------------------

.. toctree::
    :maxdepth: 2

    cases/regularBox/jupyter/VisualizeResults.ipynb


