.. include:: /include.rst_

Digital Rock Physics - Effective Permeability of Rocks
======================================================

Theoretical background
-------------------------------------

Flow in porous media occurs on different scales and differing theoretical frameworks are used to resolve flow on those scales. 
**Direct simulations** are based on solving the Navier-Stokes equation using a discretized representation(mesh) of the pore space. These simulations provide a complete description of how fluids flow through the pore network. Naturally, such simulations are computationally demanding as the nano- to micro-meter sized pores need to be resolved. In its simples incrompressible and stead-state form, it can be written as:

.. math::
    :label: eq:continuity 

    \nabla \cdot \vec{U} = 0
    
.. math::
    :label: eq:mom_con
    
    \nabla \cdot (\vec{U} \vec{U}) = \nabla \cdot (\nu \nabla \vec{U}) - \nabla p

Where :eq:`eq:continuity` is the continuity equation (based on mass conservation) and :eq:`eq:mom_con` is the momentum balance equation. :math:`\nu` is the kinematic viscosity and :math:`p` the kinematic pressure.

.. tip::
    There are many books and resources out there on the derivation of the Navier-Stokes equations, including many in an OpenFOAM context. We really like `Tobias Holtzmann's book <https://holzmann-cfd.com/community/publications/mathematics-numerics-derivations-and-openfoam>`_. Have a look! There is also the new `book by CFD direct <https://doc.cfd.direct/notes/cfd-general-principles/contents>`_.


An alternative are **Continuum simulation**, which effective properties (permeability) are used for representative control volumes are prescribed and *Darcy's law* is used to compute the flow.

.. math::
    :label: eq:darcy 

    \vec{U} = - \frac{k}{\mu_f} (\nabla p -\rho \vec{g})

Here :math:`\vec{U}` is the so-called Darcy velocity, which is the effective flow through a porous media (we will get to this later but this is actual velocity weighted by the porosity as flow only occurs in the pores). The involved constants are viscosity, :math:`\mu_f`, and permeability, :math:`\k` . Extensive theoretical work helped to elucidate the relationship between Navier-Stokes flow and Darcy flow, which was first derived as a phenomelogical equation.

Note that **direct** and **continuum** simulatons are not *per se* related to absolute scales (like mm vs meters) but simply refer to differing relative length scales / processes that are resolved.

.. tip::
    There is a great pdf on the relationship between Navier-Stokes and Darcy flow on the webpages of `Cyprien Soulaine <https://www.cypriensoulaine.com/publications-copy>`_.


.. admonition:: Further info

    Have a look at the `digital rock portal <https://www.digitalrocksportal.org>`_ to get a better idea of what Digital Rock Physics is all about! 
