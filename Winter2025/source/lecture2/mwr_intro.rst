Introduction: Method of Weighted Residuals (MWR)
=================================================

.. admonition:: Credit!

    This chapter is based on lecture notes originally prepared by Prof. Keith A. Woodbury of The University of Alabama. Further details can also be found in the book by :cite:`grandin1991`.

Suppose we have a linear differential operator :math:`D` acting on a function :math:`u`
to produce a function :math:`p`. 

.. math::
    :label: eq:mwr_basic

    D(u(x))=p(x). 

We wish to approximate :math:`u` by a functions :math:`\tilde{u}`, which is a linear
combination of basis functions chosen from a linearly independent set. That
is, 

.. math::
    :label: eq:mwr_basic_2
    
    u\cong \tilde{u}=\sum_{i=1}^na_i\varphi _i

where :math:`a` are the coefficients of the :math:`n` basis functions :math:`\varphi`. Now, when substituted into the differential operator, :math:`D`, the result of the operations is not, in general, :math:`p(x)`. Hence an error or residual will exist:

.. math::
    :label: eq:mwr_basic_3

    E(x)=R(x)=D(\tilde{u}(x))-p(x)\neq 0. 

The notion in the MWR is to force the residual to zero in some average sense
over the domain. That is

.. math::
    :label: eq:mwr_basic_4

    \int_XR(x)W_idx=0\ \ \ \ \ \ \ i=1,2,...,n  \label{MWR}


where the nuber of *weight functions* :math:`W_i` is exactly equal the number
of unknow constants :math:`a_i` in :math:`\tilde{u}`. The result is a set of :math:`\ n`
algebraic equations for the unknown constants :math:`a_i`.There are several MWR sub-methods, which differ in their choices of the weighting functions :math:`W`. Here are some popular examples:

    - Collocation method.
    - Sub-domain Method.
    - Least Squares method.
    - Galerkin method.

Each of these will be explained below. Afterwards we will apply some of them in an example to solve the steady-state advection diffusion equation.

Collocation Method
-------------------

In this method, the weighting functions are taken from the family of Dirac :math:`\delta` functions in the domain. That is, :math:`W_i(x)=\delta (x-x_i)`. The Dirac :math:`\delta` function has the property that 

.. math::
    :label: eq:mwr_collocation

    \delta (x-x_i)=\left\{ 
    \begin{array}{ccc}
    1 &  & x=x_i \\ 
    0 &  & \text{otherwise}
    \end{array}
    \right. . 

Hence the integration of the weighted residual statement results in the
forcing of the residual to zero at specific points in the domain. That is,
integration of :eq:`eq:mwr_basic_4` with :math:`W_i(x)=\delta (x-x_i)` results in 

.. math::
    :label: eq:mwr_collocation_2

    R(x_i)=0 

Sub-domain Method
------------------

This method doesn't use weighting factors explicity, so it is not, strictly
speaking, a member of the Weighted Residuals family. However, it can be
considered a modification of the collocation method. The idea is to force
the weighted residual to zero not just at fixed points in the domain, but
over various subsections of the domain. To accomplish this, the weight
functions are set to unity, and the integral over the entire domain is
broken into a number of subdomains sufficient to evaluate all unknown
parameters. That is 

.. math::
    :label: eq:mwr_subdomain_1

    \int_{X}R(x)W_{i}dx=\sum_{i}\left( \int_{X_{i}}R(x)dx\right) =0\ \ \ \ \ \ \
    i=1,2,...,n 


Least Squares Method
--------------------

If the continuous summation of all the squared residuals is minimized, the
rationale behind the name can be seen. In other words, a minimum of

.. math::
    :label: eq:mwr_lsq_1

    S=\int_XR(x)R(x)dx=\int_XR^2(x)dx. 


In order to achieve a minimum of this scalar function, the derivatives of :math:`S` with respect to all the unknown parameters must be zero. That is, 

.. math::
    :label: eq:mwr_lsq_2

    \frac{\partial S}{\partial a_i} &=&0 \\
    &=&2\int_XR(x)\frac{\partial R}{\partial a_i}dx


Comparing with :eq:`eq:mwr_basic_4`, the weight functions are seen to be 

.. math::
    :label: eq:mwr_lsq_3

    W_i=2\frac{\partial R}{\partial a_i} 

However, the :math:`2` can be dropped, since it cancels out in the equation.
Therefore the weight functions for the Least Squares Method are just the
derivatives of the residual with respect to the unknown constants: 

.. math::
    :label: eq:mwr_lsq_4
    
    W_i=\frac{\partial R}{\partial a_i} 


Galerkin Method
---------------

This method may be viewed as a modification of the Least Squares Method.
Rather than using the derivative of the *residual* with respect to the
unknown :math:`a_i`, the derivative of the approximating function is used. That
is, if the function is approximated as in :eq:`eq:mwr_basic_2`, then the weight
functions are

.. math::
    :label: eq:mwr_galerkin_1

    W_i=\frac{\partial \tilde{u}}{\partial a_i} 

Note that these are then identical to the original basis functions appearing
in :eq:`eq:mwr_basic_2`

.. math::
    :label: eq:mwr_galerkin_2

    W_i=\frac{\partial \tilde{u}}{\partial a_i}=\varphi _i(x) 
