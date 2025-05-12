.. include:: /include.rst_

Digital Rock Physics - Effective Permeability of Rocks
======================================================

.. admonition:: Further info

    Have a look at the `digital rock portal <https://digitalporousmedia.org>`_ to get a better idea of what Digital Rock Physics is all about! 


Theoretical Background
-----------------------

Flow in porous media occurs at various scales, and different theoretical frameworks are used to describe flow depending on the scale. **Direct simulations** are based on solving the Navier-Stokes equations using a discretized representation (mesh) of the pore space. These simulations provide a detailed description of the flow field within the pore network. Naturally, such simulations are computationally demanding, as the nano- to micrometer-sized pores must be fully resolved. The simplified Navier-Stokes equations for incompressible, steady-state flow are:

.. math::
    :label: eq:continuity 

    \nabla \cdot \vec{U} = 0

.. math::
    :label: eq:mom_con

    \nabla \cdot (\vec{U} \vec{U}) = \nabla \cdot (\nu \nabla \vec{U}) - \nabla p

Here, :eq:`eq:continuity` is the continuity equation (expressing mass conservation), and :eq:`eq:mom_con` is the momentum balance equation. :math:`\nu` is the kinematic viscosity, and :math:`p` is the kinematic pressure.

.. tip::
    There are many resources on the derivation of the Navier–Stokes equations, including several tailored to OpenFOAM. We particularly recommend `Tobias Holzmann's book <https://holzmann-cfd.com/community/publications/mathematics-numerics-derivations-and-openfoam>`_. Also see the new `book by CFD Direct <https://doc.cfd.direct/notes/cfd-general-principles/contents>`_.
    There is also a great pdf on the relationship between Navier-Stokes and Darcy flow on the webpages of `Cyprien Soulaine <https://www.cypriensoulaine.com/publications-copy>`_.

An alternative approach is **continuum simulation**. In continuum models, volume-averaged properties are used to describe flow behavior. The two key parameters in such models are porosity and permeability. Porosity represents the void fraction within a control volume, while permeability quantifies how easily fluids can flow through it under a pressure gradient. Although these properties are correlated, the exact relationship often needs to be determined empirically or estimated using *Digital Rock Physics* methods.

Flow at the continuum scale is described by *Darcy's law*:

.. math::
    :label: eq:darcy 

    \vec{u} = - \frac{k}{\mu_f} (\nabla p - \rho \vec{g})

Here, :math:`\vec{u}` is the **Darcy velocity**, which is the volume-averaged flow velocity through the porous medium. We will discuss this further, but note that it relates to the pore velocity via :math:`\vec{u} = \epsilon \vec{v}`, where :math:`\epsilon` is the porosity and :math:`\vec{v}` the actual fluid velocity in the pores. The parameters :math:`\mu_f` and :math:`k` are the fluid's dynamic viscosity and the medium's permeability, respectively.

Extensive theoretical work has clarified the relationship between Navier-Stokes flow and Darcy flow, the latter originally proposed as a phenomenological equation.

Note that **direct** and **continuum** simulations are not *per se* related to absolute scales (like mm vs meters) but simply refer to differing relative length scales / processes that are resolved.


.. admonition:: Important note

    It is important to note that the **absolute values** of **porosity** and **permeability** in porous media are not universal constants but depend on the **choice of the Representative Elementary Volume (REV)**. Since the REV defines the scale over which microscale heterogeneities are averaged, different choices can lead to significantly different estimates of these effective properties. This scale dependence means that **porosity and permeability are not always independently observable** — observational data often reflect a combined, scale-dependent response of the medium, making it difficult to infer one property without assumptions about the other. In media with strong heterogeneities or multi-scale structures, identifying a meaningful REV is particularly challenging. These limitations should be kept in mind when interpreting or applying continuum-scale flow models like Darcy’s law.




