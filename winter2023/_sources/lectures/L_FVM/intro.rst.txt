.. include:: /include.rst_

FVM and OpenFoam
=====================================================

OpenFOAM uses the finite volume method to solve the governing partial differential equations of a problem. In the following excercise, we will derive the FV form of the Laplance equation (transient heat diffusion) and look into the details of OpenFOAM sets up the coefficient matrix, how boundary conditions are numericaly set, and how the final matrix equation is solved. 

Most of the excercise is about theory and in lecture form. We do, however, provide jupyter notebooks in the end that allows reproducing the figures and checking the output of the modified Laplace solver against the derived values. 
