.. include:: /include.rst_

.. _Advection1:

====================
The advection term
====================

In the previous sessions, we have explored in details how openFoam solves the diffusion equation. In this session, we will explore how openFoam solves the advection equation. Here is a basic advection diffusion equation:

.. math::
    \frac{\partial T}{\partial t} = \nabla \cdot \left( D \nabla T \right) - \nabla \cdot \left( \vec{U} T \right)

where :math:`T` is the transported scalar, :math:`D` is the diffusion coefficient and :math:`\vec{U}` is the velocity field. We have already discussed the diffusion term (first term onthe RHS); now it's time to discuss advection.





