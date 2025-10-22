.. _fem_formulation:

Finite Element Formulation
===========================


Weak form
---------
To cast :ref:`eqs:stokes_equation_penalty` and :ref:`eqs:stokes_equation_divergence` into a finite element formulation, we need to define the trial and test spaces.  We follow the formulation in :cite:`Dabrowski2008` and use the the seven-node Crouzeix-Raviart triangle with quadratic velocity shape functions enhanced by a cubic bubble function and discontinuous linear interpolation for the pressure field. This element is stable and the velocity and pressure approximations are independent, which leads to the so-called mixed formulation of the finite element method.

.. math::
   :label: eq:shape_function_stokes

   \tilde{u}_x(x,y) = \sum_{i=1}^{\text{nnod}} N_i(x,y) u_x^i \\
   \tilde{u}_y(x,y) = \sum_{i=1}^{\text{nnod}} N_i(x,y) u_y^i \\
   \tilde{p}(x,y) = \sum_{i=1}^{\text{np}} \Pi_i(x,y) p^i


We can now use the Galerkin method to write the governing equations in their FEM weak form. These were the equations as derived in the previous sections:


.. math::
    :label: eqs:stokes_equation_penalty

    \begin{split}
    \frac{\partial}{\partial x} \left ( \eta \left ( \frac{4}{3} \frac{\partial u_x}{\partial x} - \frac{2}{3}\frac{\partial u_y}{\partial y} \right ) \right ) + \frac{\partial}{\partial y} \left ( \eta \left ( \frac{\partial u_x}{\partial y} + \frac{\partial u_y}{\partial x} \right ) \right ) - \frac{\partial P}{\partial x} + \rho g_x = 0\\
    \frac{\partial}{\partial y} \left ( \eta \left ( \frac{4}{3} \frac{\partial u_y}{\partial y} - \frac{2}{3}\frac{\partial u_x}{\partial x} \right ) \right ) + \frac{\partial}{\partial x} \left ( \eta \left ( \frac{\partial u_x}{\partial y} + \frac{\partial u_y}{\partial x} \right ) \right ) - \frac{\partial P}{\partial y} + \rho g_y = 0\\\
    \end{split}


And the continuity equation:

.. math::
    :label: eqs:stokes_equation_penalty_divergence

    \frac{\partial u_x}{\partial x} + \frac{\partial u_y}{\partial y} = -\frac{p}{\kappa}

We now use the shapes functions :ref:`eq:shape_function_stokes` to write the weak form of the equations. The weak form of the equations is obtained by multiplying the equations by a test function and integrating over the domain. The again use the Galerkin method, so that test functions are the same as the shape functions. We also do an integration by parts on the velocity equations. The weak form of the equations is then given by:

.. math::
   :label: eq:weak_form

   \int_{\Omega} \left[ \frac{\partial N_i}{\partial x} \left( \eta \left( \frac{4}{3} \frac{\partial (N_j u_{x,j})}{\partial x} - \frac{2}{3} \frac{\partial (N_j u_{y,j})}{\partial y} \right) \right) + \frac{\partial N_i}{\partial y} \left( \eta \left( \frac{\partial (N_j u_{x,j})}{\partial y} + \frac{\partial (N_j u_{y,j})}{\partial x} \right) \right) - \frac{\partial N_i}{\partial x} (\Pi_j p_j) - N_i f_x \right] \, d\Omega = 0 \\

   \int_{\Omega} \left[ \frac{\partial N_i}{\partial y} \left( \eta \left( \frac{4}{3} \frac{\partial (N_j u_{y,j})}{\partial y} - \frac{2}{3} \frac{\partial (N_j u_{x,j})}{\partial x} \right) \right) + \frac{\partial N_i}{\partial x} \left( \eta \left( \frac{\partial (N_j u_{x,j})}{\partial y} + \frac{\partial (N_j u_{y,j})}{\partial x} \right) \right) - \frac{\partial N_i}{\partial y} (\Pi_j p_j) - N_i f_y \right] \, d\Omega = 0 \\

   \int_{\Omega} \Pi_i \left( \frac{\partial (N_j u_{x,j})}{\partial x} + \frac{\partial (N_j u_{y,j})}{\partial y} \right) \, d\Omega = 0


Note that we obmitted the boundary terms and flipped the sign again after the integration by parts for better readability. The boundary terms contain the stresses acting on the domain. If we have stress free conditions, the terms are zero. If we want to apply a stress, we need to add the stress to the boundary terms and remember the sign flip.


Element matrices
---------------------------------

The B matrix
^^^^^^^^^^^^

We again split the integral in a sum of integrals over our triangular elements. We can then start to write the weak form in matrix form. We start with the matrix :math:`B` that turns velocities into strain rates:

.. math::
   :label: eq:B_matrix

   \begin{pmatrix}
   \dot{\epsilon}_{xx}(x,y) \\
   \dot{\epsilon}_{yy}(x,y) \\
   \dot{\epsilon}_{xy}(x,y)
   \end{pmatrix}
   = \mathbf{B}(x,y) u^e = \begin{pmatrix}
   \frac{\partial N_1}{\partial x}(x,y) & 0 & \cdots \\
   0 & \frac{\partial N_1}{\partial y}(x,y) & \cdots \\
   \frac{\partial N_1}{\partial y}(x,y) & \frac{\partial N_1}{\partial x}(x,y) & \cdots
   \end{pmatrix}
   \begin{pmatrix}
   u_x^1 \\
   u_y^1 \\
   \vdots
   \end{pmatrix}

Note how the nodal velocities are ordered in the matrix. 

The D matrix
^^^^^^^^^^^^

We then use the matrix :math:`D` to turn strain rates into stresses. Remember that deviatoric stresses are are related to deviatoric strain rates by the viscosity:

.. math::
   :label: eq:deviatoric_strain_rates

   \begin{aligned}
   \epsilon_{xx}' &= \epsilon_{xx} - \frac{1}{3}(\epsilon_{xx} + \epsilon_{yy} + \epsilon_{zz}) \\
   \epsilon_{yy}' &= \epsilon_{yy} - \frac{1}{3}(\epsilon_{xx} + \epsilon_{yy} + \epsilon_{zz}) \\
   \epsilon_{xy}' &= \epsilon_{xy}
   \end{aligned}

For a 2D plane strain condition, :math:`\epsilon_{zz} = 0`:

.. math::

   \begin{aligned}
   \epsilon_{xx}' &= \epsilon_{xx} - \frac{1}{3}(\epsilon_{xx} + \epsilon_{yy}) \\
   \epsilon_{yy}' &= \epsilon_{yy} - \frac{1}{3}(\epsilon_{xx} + \epsilon_{yy}) \\
   \epsilon_{xy}' &= \epsilon_{xy}
   \end{aligned}

Deviatoric strain rates can now be written as:

.. math::

   \begin{pmatrix}
   \epsilon_{xx}' \\
   \epsilon_{yy}' \\
   \epsilon_{xy}'
   \end{pmatrix}
   =
   \begin{pmatrix}
   \frac{2}{3} & -\frac{1}{3} & 0 \\
   -\frac{1}{3} & \frac{2}{3} & 0 \\
   0 & 0 & 1
   \end{pmatrix}
   \begin{pmatrix}
   \epsilon_{xx} \\
   \epsilon_{yy} \\
   \epsilon_{xy}
   \end{pmatrix}


Now we can use viscosity to turn deviatoric strain rates into deviatoric stresses:

.. math::

   \begin{pmatrix}
   \sigma_{xx}' \\
   \sigma_{yy}' \\
   \sigma_{xy}'
   \end{pmatrix}
   =
   2\eta
   \begin{pmatrix}
   \frac{2}{3} & -\frac{1}{3} & 0 \\
   -\frac{1}{3} & \frac{2}{3} & 0 \\
   0 & 0 & 1
   \end{pmatrix}
   \begin{pmatrix}
   \epsilon_{xx} \\
   \epsilon_{yy} \\
   \epsilon_{xy}
   \end{pmatrix}


And finally get the D matrix:

.. math::
   :label: eq:D_matrix

   \mathbf{D} =
   \eta \begin{pmatrix}
   \frac{4}{3} & -\frac{2}{3} & 0 \\
   -\frac{2}{3} & \frac{4}{3} & 0 \\
   0 & 0 & 1
   \end{pmatrix}


The coefficient/stiffness matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can now write the first two terms of the weak form in matrix form:

.. math::
   :label: eq:integral_B_TDBu

   A = \int_{\Omega_e} \eta \mathbf{B}^T \mathbf{D} \mathbf{B} \, d\Omega_e


The third term of the weak form can be written as:

.. math::
   :label: eq:integral_B_TDBu

   Q^T = -\int_{\Omega_e} \mathbf{B_{vol}}^T \Pi^T \, d\Omega_e


We can now write the final element matrix as:

.. math::
   :label: eq:stiffness_matrix

   \mathbf{K}^e = \begin{pmatrix}
   \mathbf{A} & \mathbf{Q}^T \\
   \mathbf{Q} & -\kappa^{-1} \mathbf{M}
   \end{pmatrix}
   = \int_{\Omega^e} \begin{pmatrix}
   \mu^e \mathbf{B}^T \mathbf{D} \mathbf{B} & -\mathbf{B}_{\text{vol}}^T \Pi^T \\
   -\Pi \mathbf{B}_{\text{vol}} & -\kappa^{-1} \Pi \Pi^T
   \end{pmatrix} \, dx \, dy

Note that we flipped the sign again on the pressure equation to make the coefficient matrix symmetric.


Static condensation
-------------------

Since we use a discontinuous pressure formation, we can then eliminate the pressure field from the equations. This is called static condensation. We can write the pressure field in terms of the velocity field:

.. math::
   :label: eq:pressure_static_condensation

   p(x, y) = \kappa \Pi^T(x, y) \mathbf{M}^{-1} \mathbf{Q} \mathbf{u}^e


We can now substitute the pressure field in the velocity equations to get the final coefficient matrix:

.. math::
   :label: eq:matrix_equation

   \left( \mathbf{A} + \kappa \mathbf{Q}^T \mathbf{M}^{-1} \mathbf{Q} \right) \mathbf{u}^e = \mathbf{f}^e


This is the final matrix equation that we need to solve for the velocity field.


