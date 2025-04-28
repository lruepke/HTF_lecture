.. include:: /include.rst_

Digital Rock Physics - Effective Permeability of Rocks
======================================================

.. admonition:: Further info

    Have a look at the `digital rock portal <https://www.digitalrocksportal.org>`_ to get a better idea of what Digital Rock Physics is all about! 


Theoretical background
-------------------------------------

Flow in porous media occurs on different scales and differing theoretical frameworks are used to resolve flow on those scales. **Direct simulations** are based on solving the Navier-Stokes equation using a discretized representation (mesh) of the pore space. These simulations provide a complete description of the flow field within the pore network. Naturally, such simulations are computationally demanding as the nano- to micrometer-sized pores need to be resolved. The simplified Navier-Stokes equations for incompressible steady-state flow can be written as:

.. math::
    :label: eq:continuity 

    \nabla \cdot \vec{U} = 0
    
.. math::
    :label: eq:mom_con
    
    \nabla \cdot (\vec{U} \vec{U}) = \nabla \cdot (\nu \nabla \vec{U}) - \nabla p

Where :eq:`eq:continuity` is the continuity equation (based on mass conservation) and :eq:`eq:mom_con` is the momentum balance equation. :math:`\nu` is the kinematic viscosity and :math:`p` the kinematic pressure.

.. tip::
    There are many books and resources out there on the derivation of the Navier-Stokes equations, including many in an OpenFOAM context. We really like `Tobias Holtzmann's book <https://holzmann-cfd.com/community/publications/mathematics-numerics-derivations-and-openfoam>`_. Have a look! There is also the new `book by CFD direct <https://doc.cfd.direct/notes/cfd-general-principles/contents>`_.


An alternative are **Continuum simulation**. In continuum simulations, averaged properties of control volumes are used to describe flow. The two key parameters of such continuum simulations are the porosity and the permeability of the control volumes. Porosity is the void space in each control volume and permeability describes how well fluids are transported along a pressure gradient. Both properties are  correlated but the exact relationship often needs to be emperically determined, or predicted by *Digital Rock Physics* methods. Flow on the continuum-scale is described by *Darcy's law*.

.. math::
    :label: eq:darcy 

    \vec{u} = - \frac{k}{\mu_f} (\nabla p -\rho \vec{g})

Here :math:`\vec{u}` is the so-called Darcy velocity, which is the average flow through a porous media (we will get to this later but this is pore velocity weighted by the porosity :math:`\vec{u}=\epsilon \vec{v}` , with :math:`\epsilon` being the porosity and :math:`\vec{v}` the fluid's pore velocity). The involved constants are viscosity, :math:`\mu_f`, and permeability, :math:`k` . Extensive theoretical work helped to elucidate the relationship between Navier-Stokes flow and Darcy flow, which was first derived as a phenomelogical equation.

.. tip::
    There is a great pdf on the relationship between Navier-Stokes and Darcy flow on the webpages of `Cyprien Soulaine <https://www.cypriensoulaine.com/publications-copy>`_.

Note that **direct** and **continuum** simulations are not *per se* related to absolute scales (like mm vs meters) but simply refer to differing relative length scales / processes that are resolved.



