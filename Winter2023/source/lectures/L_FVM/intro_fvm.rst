* STart from governing equations
* ask how we can solve them
* explain that openFoam uses finite volume method, and that we need to discretize the equations on the computational mesh.
* Explain finite volume method: 

The Finite Volume Method (FVM) is a popular numerical technique used in computational fluid dynamics (CFD) and other areas of computational physics to solve partial differential equations (PDEs) that describe physical phenomena. It's particularly well-suited for problems involving fluid flow, heat transfer, and associated physical processes. Here's an overview of how FVM works:

Basic Concept
Domain Discretization: The computational domain (the area or volume where the physical problem occurs) is divided into small, discrete control volumes (cells). These control volumes cover the entire domain without overlapping.
Integral Form of PDEs: FVM operates on the integral form of the governing PDEs, as opposed to the differential form used in methods like the Finite Difference Method (FDM). The equations are integrated over each control volume.
Flux Calculations: The key aspect of FVM is the calculation of fluxes of the conserved quantities (like mass, momentum, energy) across the faces of each control volume. These fluxes are used to determine how these quantities change within each control volume over time.
Boundary Conditions: The method also incorporates boundary conditions, which define how the field behaves at the boundaries of the computational domain (like the walls of a pipe, the surface of an object, etc.).
Steps in FVM
Represent the Physical Quantities: Variables like velocity, pressure, temperature are represented at discrete locations in each control volume - either at the cell centers or at the cell faces.
Discretize the Governing Equations: The integral form of the governing PDEs is discretized for each control volume. This involves expressing the rate of change of a quantity in a control volume in terms of the fluxes through its faces.
Approximate the Fluxes: The fluxes through the control volume faces are approximated using values of the variables at the centers or faces of the cells. Various schemes (like upwind, central difference) can be used for this approximation.
Solve the Discretized Equations: The discretized equations form a system of algebraic equations. These equations are solved iteratively to find the field values (like pressure, velocity, etc.) throughout the domain.
Iterate to Steady-State or March in Time: For steady-state problems, iterations continue until the solution converges. For transient problems, the solution is marched forward in time, updating the variables at each time step.
Advantages
Conservation: FVM ensures conservation of physical quantities locally (in each control volume) and globally (across the entire domain).
Flexibility: It can handle complex geometries and irregular meshes.
Applicability: Suitable for a wide range of problems in fluid dynamics and related fields.
Applications
Fluid Flow: Modeling airflow around objects, water flow in channels, etc.
Heat Transfer: Simulating heat distribution and transfer in various media.
Reactive Flows and Combustion: Used in chemical process simulations and combustion engines.
In summary, FVM is a versatile and powerful tool for solving complex physical problems numerically. It is particularly notable for its ability to handle irregular geometries and ensure the local and global conservation of physical quantities.



Openfoam uses what is called a cell centered finite volume method. This means that the variables are stored at the center of the cell. The fluxes are calculated at the cell faces. The fluxes are then used to update the variables at the cell centers.

Let's look at an example problem that solves for temperature diffusion, using the laplacianFoam solver. The governing equation is the heat diffusion equation, which is given by:

Let's start by making the mesh.




