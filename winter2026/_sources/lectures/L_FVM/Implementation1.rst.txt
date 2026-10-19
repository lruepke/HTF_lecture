.. include:: /include.rst_

.. _Implementation1:

.. tip::

   Before we start the exercise 1, please download (:download:`test_laplacianFoam <cases/test_laplacianFoam.zip>`) the modified Laplacian solver first.
   Then put it in the shared folder and compile (:code:`wmake`) it in the docker container.
   This modified solver writes out detailed information on matrix coefficients into the log file, so that it is possible to check the derived values against the computed values. 

==============
OpenFoam implementation
==============

Ok, now let's do some practical tings, 
(0) **generate mesh**; 
(1) **read mesh and do some useful calculation**, e.g. cell volume, face area, ..., etc; 
(2) **discretize Laplacian term and get coefficients**; 
(3) **discretize transient term and get additional coefficients**; 
(4) **construct the final coefficient matrix and RHS**; 
(5) **solve the system of algebraic equations**; 
(6) **write solution to file**.

.. admonition:: Goal

   Deeply look into OpenFOAM and understand how it works!

In order to better understand OpenFOAM's logic and its work flow, we have to look at the basic structure of the `source code <https://cpp.openfoam.org/v6/laplacianFoam_8C_source.html>`_ of a basic solver, :code:`laplacianFoam`.

.. tab:: laplacianFoam.C

   .. code-block::  cpp
      :linenos: 
      :emphasize-lines: 7, 9, 17, 19, 21, 23
      :caption: Source code of laplacianFoam
      :name: lst:source_laplacianFoam

      #include "fvCFD.H"               // Basic head file of OF
      #include "simpleControl.H"       // Basic head file of OF
      int main(int argc, char *argv[]) // Typical c++ main control function
      {
         #include "setRootCaseLists.H" // Do some case/file path-related thing
         #include "createTime.H"       // Create a time object: read controlDict, ...
         #include "createMesh.H"       // (2) Create mesh object: read mesh and do some useful calculation
         simpleControl simple(mesh);   // Time loop control object
         #include "createFields.H"     // (3) Read input data: T and D in the PDE 

         Info<< "\nCalculating temperature distribution\n" << endl;
         while (simple.loop(runTime)) // do time loop
         {
            Info<< "Time = " << runTime.timeName() << nl << endl;
            while (simple.correctNonOrthogonal())  // Non-orthogonal correction loop for the unstructured/non-orthogonal mesh
            {
               fvScalarMatrix TEqn                 // (5) construct the final coefficient matrix, include RHS (.source)
                  (
                     fvm::ddt(T)                   // (4) Discretization of the transient term, return a fvMatrix object
                     - 
                     fvm::laplacian(DT, T)         // (3) Discretization of the Laplacian, return a fvMatrix object
                  );
                  TEqn.solve();                    // (6) Solve the system of algebraic equations
            }
            #include "write.H"                     // (7) Save solution to file.
            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                  << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                  << nl << endl;
         }
         Info<< "End\n" << endl;
         return 0;
      }

   .. list-table:: Comparison between equation and code
      :header-rows: 0

      *  - Term
         - Equation
         - Code
         - Coefficients(ldu:b)
         - Reference
      *  - Transient
         - :math:`\frac{\partial T}{\partial t}`
         - :code:`fvm::ddt(T)`
         - d(:math:`FluxC`), b(:math:`FluxC^o`)
         - :eq:`eq:fvm_firstOrder_Euler_exp_coef`
      *  - Laplacian
         - :math:`\nabla \cdot D \nabla T`
         - :code:`fvm::laplacian(DT, T)`
         - d(:math:`\sum\limits_{F\sim NB(C)}a_F`), u(:math:`a_F` only internal faces)
         - :eq:`eq:fvm_matrix_form_internalCell`, :eq:`eq:fvm_laplacian_coeff_boundary_fixedvalue`, :eq:`eq:fvm_laplacian_coeff_boundary_fixedflux`

.. tab:: createFields.H

   .. code-block::  cpp
      :linenos:
      :emphasize-lines: 0
      :caption: createFields.H
      :name: lst:createFields

      volScalarField T //Read input data of T, include ICs and BCs
      (
         IOobject
         (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
         ),
         mesh
      );

      IOdictionary transportProperties // Read dictionary file of transportProperties
      (
         IOobject
         (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
         )
      );
      dimensionedScalar DT // Read the diffusivity D (constant) from transportProperties object
      (
         transportProperties.lookup("DT")
      );

.. tab:: Input data 1

   .. code-block::  foam
      :linenos:
      :emphasize-lines: 0
      :caption: constant/transportProperties
      :name: lst:transportProperties

      FoamFile
      {
         version     2.0;
         format      ascii;
         class       dictionary;
         location    "constant";
         object      transportProperties;
      }
      DT              DT [0 2 -1 0 0 0 0] 4e-05;

.. tab:: Input data 2

   .. code-block::  foam
      :linenos:
      :emphasize-lines: 0
      :caption: 0/T
      :name: lst:0_T

      FoamFile
      {
         version     2.0;
         format      ascii;
         class       volScalarField;
         object      T;
      }
      dimensions      [0 0 0 1 0 0 0];
      internalField   uniform 273;
      boundaryField
      {
         left
         {
            type            fixedValue;
            value           uniform 273;
         }
         right
         {
            type            fixedValue;
            value           uniform 573;
         }
         "(top|bottom)"
         {
            type            zeroGradient;
         }
      }


Step0, Generate mesh
---------------------------

.. tab:: Generate mesh

   Nothing special, just a meshing process. 
   Here we use the OpenFOAM utility :code:`blockMesh` to generate a regular mesh (:download:`Regular box case <cases/regularBox.zip>`) same as :numref:`fig:polyMesh_regularBox`.

.. tab:: OpenFOAM script

   .. code-block::  foam
      :linenos:
      :emphasize-lines: 0
      :caption: blockMeshDict of a regular box mesh we shown above
      :name: lst:blockMesh_regularBox

      /*--------------------------------*- C++ -*----------------------------------*\
      | =========                 |                                                 |
      | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
      |  \\    /   O peration     | Version:  5                                     |
      |   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
      |    \\/     M anipulation  |                                                 |
      \*---------------------------------------------------------------------------*/
      FoamFile
      {
         version     2.0;
         format      ascii;
         class       dictionary;
         object      blockMeshDict;
      }
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

      convertToMeters 1;
      xmin -0.05;    //variable definition
      xmax 0.05;
      ymin -0.015;
      ymax 0.015;
      Lz 0.01;
      vertices    //vertices definition
      (
         ($xmin      $ymin   0)  //coordinate of vertex 0
         ($xmax    $ymin   0)  //coordinate of vertex 1
         ($xmax    $ymax   0)  //coordinate of vertex 2
         ($xmin      $ymax   0)  //coordinate of vertex 3
         ($xmin      $ymin   $Lz)//coordinate of vertex 4
         ($xmax    $ymin   $Lz)//coordinate of vertex 5
         ($xmax    $ymax   $Lz)//coordinate of vertex 6
         ($xmin      $ymax   $Lz)//coordinate of vertex 7
      );
      blocks
      (
         hex (0 1 2 3 4 5 6 7) (10 3 1) simpleGrading (1 1 1)
      );
      boundary
      (
         left    //patch name
         {
            type patch ;
            faces   //face list
            (
                  (0 4 7 3)
            );
         }
         right
         {
            type patch;
            faces
            (
                  (2 6 5 1)
            );
         }
         top
         {
            type patch;
            faces
            (
                  (3 7 6 2)
            );
         }
         bottom
         {
            type patch;
            faces
            (
                  (1 5 4 0)
            );
         }
         frontAndBack    //patch name
         {
            type empty; 
            faces //face list
            (
                  (0 3 2 1)   //back face
                  (4 5 6 7)   //front face
            );
         }
      );
      // ************************************************************************* //

Ok, let's exploring the main steps of the :code:`laplacianFoam`.
Please download (:download:`test_laplacianFoam <cases/test_laplacianFoam.zip>`) and do the following practice/debug steps in :code:`test_laplacianFoam.C`.

Step 1, Read mesh and input field
-------------------------------------

.. tab:: Read data

   The basic properties of cells and faces, e.g. area, volume, face normal vector, will be evaluated after mesh reading, all these processes are happened in the mesh object. 
   It means that after calling :code:`creatFields.H` all these properties are evaluated and stored in the :code:`mesh` object. Of course the temperature field object :code:`T` (volScalarField) with BCs and ICs, and thermal diffusivity :code:`DT` are also initialized from input data after calling :code:`creatFields.H`.
   **It should be noted** that the part of Laplacian discretization coefficients are also calculated in this step. If the mesh is not changed during simulation time, the mesh related coefficients just need to be calculated one time.

   .. list-table:: Mesh- and Field-related coefficients
      :header-rows: 0

      *  - Object
         - Equation
         - Code
         - Faces
         - Reference
      *  - :code:`mesh`
         - :math:`1/\delta_{C\leftrightarrow F}`
         - :code:`mesh.deltaCoeffs()`
         - All faces
         - :eq:`eq:fvm_matrix_form_internalCell`, :eq:`eq:fvm_laplacian_coeff_boundary_fixedvalue`, :eq:`eq:fvm_laplacian_coeff_boundary_fixedflux`
      *  - :code:`T`
         - Dependents on BCs type
         - :code:`gradientInternalCoeffs` (diagonal), :code:`gradientBoundaryCoeffs` (source)
         - Boundary patch/faces
         - :eq:`eq:fvm_laplacian_coeff_boundary_fixedvalue`, :eq:`eq:fvm_laplacian_coeff_boundary_fixedflux`

.. tab:: Access mesh and field properties

   .. code-block::  cpp
      :linenos:
      :emphasize-lines: 0
      :caption: Access mesh properties and field boundary properties
      :name: lst:access_mesh

      // (1). read data
      #include "createMesh.H"
      #include "createFields.H"
      // 1.1 access internal faces
      Info<<"\n\nAccess internal mesh"<<endl;
      surfaceVectorField Cf = mesh.Cf();
      surfaceVectorField Sf = mesh.Sf();
      surfaceScalarField S = mesh.magSf();
      forAll(Cf, iface)
      {
         Info<<iface<<": face center "<<Cf[iface]<<endl;
         Info<<iface<<": face area vector "<<Sf[iface]<<endl;
         Info<<iface<<": face area "<<S[iface]<<endl;
         Info<<iface<<": face delta coeff "<<mesh.deltaCoeffs()[iface]<<endl;
         Info<<iface<<": coeff(D*magSf*deltacoeff) "<<mesh.deltaCoeffs()[iface]*DT*S[iface]<<endl;
      }
      // 1.2. access boundary mesh
      Info<<"\n\nAccess boundary mesh"<<endl;
      const fvBoundaryMesh& boundaryMesh = mesh.boundary(); 
      forAll(boundaryMesh, patchI)
      {
         const fvPatch& patch = boundaryMesh[patchI];
         forAll(patch, faceI)
         {
            Info<<"Patch "<<patch.name()<<" face "<<faceI<<": face center "<<patch.Cf()[faceI]<<endl;
            Info<<"Patch "<<patch.name()<<" face "<<faceI<<": face area vector "<<patch.Sf()[faceI]<<endl;
            Info<<"Patch "<<patch.name()<<" face "<<faceI<<": face area "<<patch.magSf()[faceI]<<endl;
            Info<<"Patch "<<patch.name()<<" face "<<faceI<<": face delta coeff "<<patch.deltaCoeffs()[faceI]<<endl;
            Info<<"Patch "<<patch.name()<<" face "<<faceI<<": owner cell "<<patch.faceCells()[faceI]<<endl;
         } 
      }
      // 1.3. access boundary field, boundary field coefficients, 
      forAll(T.boundaryField(), patchI)
      {
         Info<<"Boundary patch: "<<mesh.boundary()[patchI].name()<<endl;
         Info<<"Is coupled ? "<<mesh.boundary()[patchI].coupled()<<endl;
         Info<<"gradientInternalCoeffs of field T "<<endl;
         Info<<T.boundaryField()[patchI].gradientInternalCoeffs()<<endl; //Diagonal coeff [A]
         Info<<"gradientBoundaryCoeffs of field T "<<endl;
         Info<<T.boundaryField()[patchI].gradientBoundaryCoeffs()<<"\n"<<endl; //source coeff, [B]
      }

.. tab:: Internal mesh

   .. literalinclude::  /_static/log_case/L_FVM/regularMesh_internal.txt
      :language: bash
      :linenos:
      :lines: 0-
      :emphasize-lines: 0
      :caption: Internal mesh properties of the regular mesh shown in :numref:`fig:polyMesh_regularBox`.
      :name: lst:log_internalMesh

.. tab:: Boundary mesh

   .. literalinclude::  /_static/log_case/L_FVM/regularMesh_boundary.txt
      :language: bash
      :linenos:
      :lines: 0-
      :emphasize-lines: 0
      :caption: Boundary mesh properties of the regular mesh shown in :numref:`fig:polyMesh_regularBox`.
      :name: lst:log_boundaryMesh

.. tab:: Boundary T

   .. literalinclude::  /_static/log_case/L_FVM/regularMesh_T_boundary.txt
      :language: bash
      :linenos:
      :lines: 0-
      :emphasize-lines: 0
      :caption: Boundary properties of field T of the regular mesh shown in :numref:`fig:polyMesh_regularBox`.
      :name: lst:log_boundaryT

.. tab:: Internal cell
   :new-set:

   .. figure:: /_figures/Coordinate_delta_internalcell_regularBox.*
      :align: center
      :width: 100 %
      :name: fig:deltaCoeff_InternalCell

      Information of internal cell (:math:`C_{12}`)

.. tab:: Boundary cell

   .. figure:: /_figures/Coordinate_delta_boundary_regularBox.*
      :align: center
      :width: 100 %
      :name: fig:deltaCoeff_BoundaryCell

      Information of boundary cell (:math:`C_{19}`)

.. _OF_fvmLaplacian:

Step 2, discretize Laplacian term
-------------------------------------

.. tab:: discretize Laplacian term

   Because discretization coefficients matrix of Laplacian term is a symmetric matrix, so :code:`fvm::Laplacian(DT, T)` will return a fvMatrix object only has diagonal and upper. 
   What :code:`fvm::Laplacian` actually did is evaluate (1) :math:`a_F` (see :eq:`eq:fvm_matrix_form_internalCell`) for each internal faces, (2) :math:`a_C` for each cells, which is the negative summation of :math:`a_F`, (3) store the field BCs-related coefficients as :code:`internalCoeffs` and :code:`boundaryCoeffs`, respectively.
   All of these are implemented in gaussLaplacianScheme.C_ . **Note that** the :code:`Gaussian` scheme is the only choice for Laplacian discretization in OF. 

.. tab:: Access Laplacian coeffs 

   .. code-block::  cpp
      :linenos:
      :emphasize-lines: 0
      :caption: Access fvm::Laplacian discretization
      :name: lst:access_Laplacian

      fvScalarMatrix Laplacian(fvm::laplacian(DT, T));
      Info<<"fvm::laplacian(DT, T): "<<"\n"
         <<"\tLower"<<Laplacian.lower()<<"\n"
         <<"\tDiagonal"<<Laplacian.diag()<<"\n"
         <<"\tUpper"<<Laplacian.upper()<<"\n"
         <<"\tinternalCoeffs"<<Laplacian.internalCoeffs()<<"\n"
         <<"\tboundaryCoeffs"<<Laplacian.boundaryCoeffs()<<"\n"
         <<"\tSource"<<Laplacian.source()<<"\n"
         <<endl;

   .. tip::

      The source coefficients come from BCs are not stored in the :code:`.source()`, but in :code:`boundaryCoeffs`.
      So if you print :code:`Laplacian.source()`, it will display zero.
      The BCs-related source will be added into the :code:`TEqn.source` when :code:`TEqn.solve()` is calling. 
      There is a protected member function named :code:`addBoundarySource` will be called in :code:`solve()` function.

.. tab:: Source code of Laplacian

   .. code-block::  cpp
      :linenos:
      :emphasize-lines: 15, 20, 39, 40
      :caption: Key source code of fvm::Laplacian (gaussLaplacianScheme.C_ before nonOrthogonal correcting)
      :name: lst:source_fvm_Laplacian

      template<class Type, class GType>
      tmp<fvMatrix<Type>>
      gaussLaplacianScheme<Type, GType>::fvmLaplacianUncorrected
      (
         const surfaceScalarField& gammaMagSf,
         const surfaceScalarField& deltaCoeffs,
         const GeometricField<Type, fvPatchField, volMesh>& vf
      )
      {
         tmp<fvMatrix<Type>> tfvm
         (
            new fvMatrix<Type>
            (
                  vf,
                  deltaCoeffs.dimensions()*gammaMagSf.dimensions()*vf.dimensions()
            )
         );
         fvMatrix<Type>& fvm = tfvm.ref();

         fvm.upper() = deltaCoeffs.primitiveField()*gammaMagSf.primitiveField();
         fvm.negSumDiag();

         forAll(vf.boundaryField(), patchi)
         {
            const fvPatchField<Type>& pvf = vf.boundaryField()[patchi];
            const fvsPatchScalarField& pGamma = gammaMagSf.boundaryField()[patchi];
            const fvsPatchScalarField& pDeltaCoeffs =
                  deltaCoeffs.boundaryField()[patchi];

            if (pvf.coupled())
            {
                  fvm.internalCoeffs()[patchi] =
                     pGamma*pvf.gradientInternalCoeffs(pDeltaCoeffs);
                  fvm.boundaryCoeffs()[patchi] =
                     -pGamma*pvf.gradientBoundaryCoeffs(pDeltaCoeffs);
            }
            else
            {
                  fvm.internalCoeffs()[patchi] = pGamma*pvf.gradientInternalCoeffs();
                  fvm.boundaryCoeffs()[patchi] = -pGamma*pvf.gradientBoundaryCoeffs();
            }
         }

         return tfvm;
      }

.. tab:: Laplacian coefficients

   .. literalinclude::  /_static/log_case/L_FVM/regularMesh_fvm_Laplacian.txt
      :language: bash
      :linenos:
      :lines: 0-
      :emphasize-lines: 58, 63, 67, 68, 77, 70
      :caption: fvm::Laplacian Coefficients  of the regular mesh shown in :numref:`fig:polyMesh_regularBox`.
      :name: lst:log_fvm_Laplacian

.. tab:: Internal cell
   :new-set:

   .. figure:: /_figures/Coordinate_Laplacian_internalcell_regularBox.*
      :align: center
      :width: 100 %
      :name: fig:LaplacianCoeff_InternalCell

      Information of internal cell (:math:`C_{12}`)

.. tab:: Boundary cell 19

   .. figure:: /_figures/Coordinate_Laplacian_boundary_C19_regularBox.*
      :align: center
      :width: 100 %
      :name: fig:LaplacianCoeff_BoundaryCell

      Information of boundary cell (:math:`C_{19}`)

.. tab:: Boundary cell 9

   .. figure:: /_figures/Coordinate_Laplacian_boundary_C9_regularBox.*
      :align: center
      :width: 100 %

      Information of boundary cell (:math:`C_{9}`)

.. tab:: Boundary cell 0

   .. figure:: /_figures/Coordinate_Laplacian_boundary_C0_regularBox.*
      :align: center
      :width: 100 %

      Information of boundary cell (:math:`C_{0}`)

.. tab:: Boundary cell 10

   .. figure:: /_figures/Coordinate_Laplacian_boundary_C10_regularBox.*
      :align: center
      :width: 100 %

      Information of boundary cell (:math:`C_{10}`)

.. tab:: Boundary cell 5

   .. figure:: /_figures/Coordinate_Laplacian_boundary_C5_regularBox.*
      :align: center
      :width: 100 %

      Information of boundary cell (:math:`C_{5}`)

.. _OF_fvmDdt:

Step 3, discretize transient term
-------------------------------------



.. tab:: discretize Laplacian term

   For implicit discretization, :code:`fvm::ddt(T)` will return a fvMatrix object contains diagonal coefficients and source.
   The coefficients depend on discretization scheme.
   For example Euler scheme, the diagonal coefficients are calculated from :eq:`eq:fvm_temporal_linear` and :eq:`eq:fvm_firstOrder_Euler_imp_coeff`.
   All these are implemented in EulerDdtScheme.C_. 

   .. tip::

      There are 8 schemes for transient discretization in OpenFOAM

      #. CoEuler
      #. CrankNicolson
      #. Euler
      #. SLTS
      #. backward
      #. bounded
      #. localEuler
      #. steadyState

.. tab:: Access fvm::ddt coeffs

   .. code-block::  cpp
      :linenos:
      :emphasize-lines: 0
      :caption: Access fvm::ddt discretization
      :name: lst:access_ddt

      Info<<"fvm::ddt(T): "<<"\n"
         <<"\tLower"<<ddt.lower()<<"\n"
         <<"\tDiagonal"<<ddt.diag()<<"\n"
         <<"\tUpper"<<ddt.upper()<<"\n"
         <<"\tinternalCoeffs"<<ddt.internalCoeffs()<<"\n" //actually this is not necessary for fvm::ddt, this is definitely equal to zero
         <<"\tboundaryCoeffs"<<ddt.boundaryCoeffs()<<"\n"
         <<"\tSource"<<ddt.source()<<"\n"
         <<endl;

.. tab:: Source code of fvm::ddt

   .. code-block::  cpp
      :linenos:
      :emphasize-lines: 19, 21, 29
      :caption: Key source code of fvm::ddt (EulerDdtScheme.C_ )
      :name: lst:source_fvm_ddt

      template<class Type>
      tmp<fvMatrix<Type>>
      EulerDdtScheme<Type>::fvmDdt
      (
         const GeometricField<Type, fvPatchField, volMesh>& vf
      )
      {
         tmp<fvMatrix<Type>> tfvm
         (
            new fvMatrix<Type>
            (
                  vf,
                  vf.dimensions()*dimVol/dimTime
            )
         );

         fvMatrix<Type>& fvm = tfvm.ref();

         scalar rDeltaT = 1.0/mesh().time().deltaTValue(); // 1/dt

         fvm.diag() = rDeltaT*mesh().Vsc(); // Vc/dt (FluxC)

         if (mesh().moving())
         {
            fvm.source() = rDeltaT*vf.oldTime().primitiveField()*mesh().Vsc0();
         }
         else
         {
            fvm.source() = rDeltaT*vf.oldTime().primitiveField()*mesh().Vsc(); // -T_old*Vc/dt
         }

         return tfvm;
      }

.. tab:: fvm::ddt coefficients

   :math:`\Delta t = 0.05\ s`

   .. literalinclude::  /_static/log_case/L_FVM/regularMesh_fvm_ddt.txt
      :language: bash
      :linenos:
      :lines: 0-
      :emphasize-lines: 0
      :caption: fvm::ddt Coefficients  of the regular mesh shown in :numref:`fig:polyMesh_regularBox`.
      :name: lst:log_fvm_ddt


Step 4, construct the final coefficient matrix and RHS
-----------------------------------------------------------

.. tab:: Final matrix

   The final coefficient matrix is constructed by simply adding the matrix of :ref:`OF_fvmLaplacian` and :ref:`OF_fvmDdt`.
   
   * The diagonal coefficients come from :math:`a_C` of Laplacian term and :math:`FluxC` of transient term
   * The off-diagonal coefficients only come from :math:`a_F (internal\ face)` of Laplacian term
   * The RHS comes from :math:`c_F` of Laplacian term when at boundary faces and :math:`FluxC^oT_C^o` of transient term.

.. tab:: Coefficients at the first time step

   .. literalinclude::  /_static/log_case/L_FVM/regularMesh_fvm_TEqn.txt
      :language: bash
      :linenos:
      :lines: 0-
      :emphasize-lines: 69, 76, 57, 66, 62, 67
      :caption: Final matrix coefficients of the regular mesh shown in :numref:`fig:polyMesh_regularBox` at the first time step.
      :name: lst:log_fvm_TEqn

Step 5, solve
-------------------------------------

For the :download:`Regular box case <cases/regularBox.zip>` case, we can use **PBiCG** solver and **DILU** preconditioner.

.. admonition:: Available preconditioner in OpenFOAM

   * **diagonal** : for symmetric & nonsymmetric matrices (not very effective)
   * **DIC** : Diagonal Incomplete Cholesky preconditioner for symmetric matrices
   * **DILU** : Diagonal Incomplete LU preconditioner for nonsymmetric matrices
   * **FDIC** : Fast Diagonal Incomplete Cholesky preconditioner
   * **GAMG** : Geometric Agglomerated algebraic MultiGrid preconditioner


.. admonition:: Available solver in OpenFOAM

   * **BICCG**: Diagonal incomplete LU preconditioned BiCG solver
   * **diagonalSolver**: Solver for symmetric and nonsymmetric matrices
   * **GAMG**: Geometric Agglomerated algebraic Multi-Grid solver
   * **ICC**: Incomplete Cholesky Conjugate Gradient solver
   * **PBiCG**: Bi-Conjugate Gradient solver with preconditioner
   * **PCG**: Conjugate Gradient solver with preconditioner
   * **smoothSolver**: Iterative solver with run-time selectable smoother

.. admonition:: Krylov-subspace solvers

   * **CG**: The Conjugate Gradient algorithm applies to systems where A is symmetric positive definite (SPD)
   * **GMRES**: The Generalized Minimal RESidual algorithm is the first method to try if A is not SPD.
   * **BiCG**: The BiConjugate Gradient algorithm applies to general linear systems, but the convergence can be quite erratic.
   * **BiCGstab**: The stabilized version of the BiConjugate Gradient algorithm.



Step 6, write
-------------------------------------




Jupyter notebook
-------------------

.. toctree::
    :maxdepth: 2

    cases/jupyter/VisualizeResults.ipynb


