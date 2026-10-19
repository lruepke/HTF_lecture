# Lecture 6: Transient and Coupled Problems - Python Implementation

## Overview

This directory contains Python implementations of the Finite Element Method (FEM) for solving transient heat transfer and coupled reaction-diffusion problems. The code has been refactored to improve maintainability, documentation, and educational clarity while fixing critical bugs.

## Recent Refactoring (2025)

### Critical Bugs Fixed

1. **Shape Function API Inconsistency** - Fixed incompatible function signatures that caused runtime errors
2. **Undefined Variable Bug** - Fixed Jacobian computation using stale derivatives in `turing_pattern.py`

### New Utility Modules

- `fem_shapes.py` - Unified shape function module with comprehensive documentation
- `fem_utils.py` - FEM assembly helpers and Jacobian computations
- `fem_integration.py` - Gauss quadrature rules for triangular elements
- `mesh_utils.py` - Mesh generation utilities
- `validation.py` - Stability checks and parameter validation

## Directory Structure

**Note:** All files are in a flat directory structure (no subdirectories). The categories below are for organizational clarity only.

```
lecture6/python/
├── README.md                          # This file
├── REFACTORING_SUMMARY.md             # Detailed refactoring documentation
├── test_imports.py                    # Import verification tests
│
├── Main Simulation Scripts:
├── turing_pattern.py                  # Turing patterns with adaptive time stepping
├── double_diff.py                     # Coupled reaction-diffusion (simpler version)
├── 2d_fem_transient_triangle.py       # Transient heat conduction
├── turing_fd_explicit.py              # Explicit FDM (for comparison)
│
├── Utility Modules (NEW - 2025):
├── fem_shapes.py                      # Shape functions for triangular elements
├── fem_utils.py                       # FEM assembly and Jacobian utilities
├── fem_integration.py                 # Gauss quadrature rules
├── mesh_utils.py                      # Mesh generation helpers
├── validation.py                      # Stability checks and validation
│
└── Deprecated Files (still work but use new modules instead):
    ├── shapes_tri.py                  # DEPRECATED → Use fem_shapes.py
    ├── shapes_tri_nnodel.py           # DEPRECATED → Use fem_shapes.py
    └── int_points_triangle.py         # DEPRECATED → Use fem_integration.py
```

## Main Scripts

**Integration Status:** ✅ All main scripts have been fully integrated with the utility modules (2025).

### 1. `turing_pattern.py`

**Purpose:** Solve coupled reaction-diffusion equations to generate Turing patterns.

**Utility Integration:**
- ✅ Uses `fem_shapes.shapes_tri()` for shape functions
- ✅ Uses `fem_integration.ip_triangle()` for Gauss quadrature
- ✅ Uses `fem_utils.compute_jacobian()` for Jacobian computation
- ✅ Uses `fem_utils.create_assembly_indices_2dof()` for assembly indices
- ✅ Uses `mesh_utils.make_box()` for mesh generation
- ✅ Uses `validation.check_cfl_condition()` for stability checks
- ✅ Uses `validation.validate_parameters()` for parameter validation

**Mathematical Model:**
```
∂A/∂t = D_A ∇²A - AB² + f(1-A)
∂B/∂t = D_B ∇²B + AB² - (k+f)B
```

**Key Features:**
- Implicit FEM with nonlinear iteration
- Adaptive time stepping
- Unstructured triangular mesh (128×128 domain)
- Gray-Scott model for pattern formation

**Parameters:**
- `D_A = 0.25` - Diffusivity of species A [m²/s]
- `D_B = 0.125` - Diffusivity of species B [m²/s]
- `f = 0.055` - Feed rate (dimensionless)
- `k = 0.062` - Kill rate (dimensionless)
- `dt = 0.25` - Time step [s]

**Output:** XDMF format for ParaView visualization

**Usage:**
```bash
python turing_pattern.py
```

### 2. `double_diff.py`

**Purpose:** Simpler version of coupled diffusion on smaller domain (5×5).

**Utility Integration:**
- ✅ Uses `fem_shapes.shapes_tri()` for shape functions
- ✅ Uses `fem_integration.ip_triangle()` for Gauss quadrature
- ✅ Uses `fem_utils.compute_jacobian()` for Jacobian computation
- ✅ Uses `fem_utils.create_assembly_indices_2dof()` for assembly indices
- ✅ Uses `mesh_utils.make_box()` for mesh generation
- ✅ Uses `validation.check_cfl_condition()` for stability checks
- ✅ Uses `validation.validate_parameters()` for parameter validation

**Differences from `turing_pattern.py`:**
- Smaller domain for faster computation
- Fixed time stepping (no adaptation)
- Good for testing and learning

### 3. `2d_fem_transient_triangle.py`

**Purpose:** Transient heat conduction in composite materials.

**Utility Integration:**
- ✅ Uses `fem_shapes.shapes_tri()` for shape functions
- ✅ Uses `fem_integration.ip_triangle()` for Gauss quadrature
- ✅ Uses `fem_utils.compute_jacobian()` for Jacobian computation
- ✅ Uses `mesh_utils.make_box()` and `make_inclusion()` for mesh generation
- ✅ Uses `validation.check_cfl_condition()` for stability checks
- ✅ Uses `validation.validate_parameters()` for parameter validation

**Mathematical Model:**
```
ρc_p ∂T/∂t = ∇·(k∇T)
```

**Key Features:**
- Multiple circular inclusions with different conductivities
- Dirichlet boundary conditions (fixed temperatures)
- Single DOF per node (temperature only)

**Parameters:**
- `k1 = 1.0` - Matrix conductivity [W/(m·K)]
- `k2 = 0.01` - Inclusion conductivity [W/(m·K)]
- `dt = 0.015` - Time step [s]

### 4. `turing_fd_explicit.py`

**Purpose:** Explicit finite difference method for Turing patterns.

**Utility Integration:**
- ✅ Uses `validation.check_cfl_condition()` for stability checks (CRITICAL for explicit methods!)
- ✅ Uses `validation.validate_parameters()` for parameter validation
- Note: This script uses finite differences (not FEM), so it doesn't need FEM utility modules

**Key Features:**
- Structured Cartesian grid (128×128)
- Periodic boundary conditions
- Very fast but limited by CFL condition
- **CFL stability check is critical** - script will warn if time step is too large

**Use Case:** Comparison with FEM for educational purposes

## Utility Modules

### `fem_shapes.py`

**Shape functions for triangular elements**

```python
from fem_shapes import shapes_tri

# Evaluate at integration point
N, dNds = shapes_tri(xi=0.333, eta=0.333, nnodel=3)

# N: shape functions, shape (nnodel,)
# dNds: derivatives w.r.t. local coords, shape (2, nnodel)
```

**Supported Elements:**
- `nnodel=3`: Linear triangles (P1)
- `nnodel=6`: Quadratic triangles (P2)

**Features:**
- Comprehensive NumPy-style docstrings
- Built-in verification functions
- Default parameter for backward compatibility

### `fem_utils.py`

**FEM assembly helpers and Jacobian computations**

```python
from fem_utils import compute_jacobian, create_assembly_indices_2dof

# Compute Jacobian for element mapping
Jac, invJ, detJ = compute_jacobian(dNds, ECOORD)

# Create assembly indices for 2-DOF system (coupled problems)
I, J = create_assembly_indices_2dof(EL2NOD, nnodel=3)
```

**Functions:**
- `compute_jacobian()` - Jacobian matrix with error checking
- `create_assembly_indices_2dof()` - Efficient index generation for coupled systems
- `assemble_sparse_matrix_2dof()` - Global matrix assembly
- `get_triangle_quadrature()` - Convenience wrapper for integration

### `fem_integration.py`

**Gauss quadrature rules for triangular elements**

```python
from fem_integration import ip_triangle

# Get 3-point quadrature (most common for linear elements)
ipx, ipw = ip_triangle(nip=3)

# ipx: integration points, shape (nip, 2)
# ipw: integration weights, shape (nip,)
```

**Supported Rules:**
- `nip=1`: 1-point (degree 1)
- `nip=3`: 3-point (degree 2) - **recommended for P1**
- `nip=6`: 6-point (degree 4)
- `nip=7`: 7-point (degree 5) - **recommended for P2**
- `nip=12, 25`: Higher-order rules
- `nip=77`: Labatto rule (includes vertices)

### `mesh_utils.py`

**Mesh generation utilities**

```python
from mesh_utils import make_box, make_inclusion
import triangle as tr

# Initialize lists
vertices, segments, regions = [], [], []

# Create outer domain
make_box(vertices, segments, regions, 0, 0, 1, 1, attribute=1)

# Add circular inclusion
make_inclusion(vertices, segments, regions,
               0.5, 0.5, 0.2, 20, attribute=100)

# Generate mesh
A = dict(vertices=vertices, segments=segments, regions=regions)
B = tr.triangulate(A, 'pq33Aa')
GCOORD = B['vertices']
EL2NOD = B['triangles']
```

**Functions:**
- `make_box()` - Rectangular domains
- `make_inclusion()` - Circular inclusions

### `validation.py`

**Stability checks and parameter validation**

```python
from validation import check_cfl_condition, validate_mesh_quality

# Check CFL stability condition
is_stable, cfl, cfl_max = check_cfl_condition(
    dt=0.01, dx=0.1, diffusivity=1e-6, method='explicit'
)

# Check mesh quality
metrics = validate_mesh_quality(GCOORD, EL2NOD)
print(f"Min angle: {metrics['min_angle']:.1f}°")
print(f"Min quality: {metrics['min_quality']:.3f}")
```

**Functions:**
- `check_cfl_condition()` - Verify time step stability
- `validate_mesh_quality()` - Check mesh quality metrics
- `check_convergence()` - Monitor iterative solver
- `validate_parameters()` - Validate input parameters

## Mathematical Background

### Weak Form of Heat Equation

The transient heat equation:
```
ρc_p ∂T/∂t = ∇·(k∇T) + Q
```

Weak form (after integration by parts):
```
∫_Ω ρc_p (∂T/∂t) φ dΩ + ∫_Ω k ∇T·∇φ dΩ = ∫_Ω Q φ dΩ + ∫_Γ q̄ φ dΓ
```

### FEM Discretization

After spatial discretization with shape functions N_i:
```
T(x,t) ≈ Σ N_i(x) T_i(t)
```

This leads to the system:
```
M dT/dt + K T = F
```

where:
- M: mass matrix
- K: stiffness matrix
- F: force vector

### Time Integration

**Implicit Euler (backward Euler):**
```
(M + Δt K) T^(n+1) = M T^n + Δt F^(n+1)
```

- Unconditionally stable
- Requires solving linear system each time step
- Used in `turing_pattern.py` and `2d_fem_transient_triangle.py`

**Explicit Euler (forward Euler):**
```
M T^(n+1) = M T^n - Δt K T^n + Δt F^n
```

- Conditionally stable (CFL condition)
- No linear solve required
- Used in `turing_fd_explicit.py`

**CFL Stability Condition:**

For explicit methods in 2D:
```
Δt ≤ (Δx)² / (4D)
```

where D is the diffusion coefficient.

## Physical Parameters

### Thermal Parameters

- **Thermal diffusivity:** α = k / (ρc_p) [m²/s]
  - Metals: ~10⁻⁴ to 10⁻⁵ m²/s
  - Rock: ~10⁻⁶ m²/s
  - Water: ~10⁻⁷ m²/s

- **Thermal conductivity:** k [W/(m·K)]
  - Copper: 400 W/(m·K)
  - Rock: 2-4 W/(m·K)
  - Air: 0.025 W/(m·K)

- **Density:** ρ [kg/m³]
- **Heat capacity:** c_p [J/(kg·K)]

### Gray-Scott Model Parameters

Parameter combinations that produce different patterns:

| Pattern | Feed rate (f) | Kill rate (k) |
|---------|---------------|---------------|
| Spots   | 0.0545        | 0.062         |
| Stripes | 0.055         | 0.062         |
| Waves   | 0.014         | 0.054         |
| Chaos   | 0.026         | 0.051         |

## Dependencies

Required Python packages:
```bash
pip install numpy scipy matplotlib triangle meshio tabulate
```

Or using the provided environment:
```bash
conda env create -f ../../environment.yaml
```

## Common Tasks

### Running a Simulation

```bash
cd /path/to/lecture6/python
python turing_pattern.py
```

Output files:
- `transient.xmf` - XDMF file for ParaView
- `transient.h5` - HDF5 data file

### Visualizing Results

**Using ParaView:**
1. Open `transient.xmf` in ParaView
2. Apply "Warp By Scalar" filter (optional)
3. Color by field variable
4. Play animation

**Using Python:**
```python
import meshio

# Read last time step
mesh = meshio.read("transient.xmf")
A = mesh.point_data['A']  # Species A concentration
B = mesh.point_data['B']  # Species B concentration
```

### Modifying Parameters

Edit the parameter section at the top of each script:

```python
# Model parameters
dt = 0.25      # Time step [s]
tottime = 1000 # Total simulation time [s]
D_A = 1/4.0    # Diffusivity A [m²/s]
D_B = 0.5/4.0  # Diffusivity B [m²/s]
f_coeff = 0.055  # Feed rate
k_coeff = 0.062  # Kill rate
```

### Creating Custom Meshes

```python
from mesh_utils import make_box, make_inclusion
import triangle as tr

vertices, segments, regions = [], [], []

# Outer domain
make_box(vertices, segments, regions, 0, 0, 2, 2, attribute=1)

# Add multiple inclusions
centers = [(0.5, 0.5), (1.5, 0.5), (1.0, 1.5)]
for cx, cy in centers:
    make_inclusion(vertices, segments, regions,
                   cx, cy, 0.2, 20, attribute=100)

# Generate mesh with quality constraints
A = dict(vertices=vertices, segments=segments, regions=regions)
B = tr.triangulate(A, 'pq30Aa0.01')  # q30: min angle 30°, a0.01: max area 0.01
```

## Troubleshooting

### Common Errors

**1. `TypeError: shapes_tri() takes 2 positional arguments but 3 were given`**

Solution: Update import to use new `fem_shapes` module:
```python
from fem_shapes import shapes_tri  # NEW
# from shapes_tri import shapes_tri  # OLD
```

**2. Simulation is unstable (oscillations or explosions)**

Check CFL condition:
```python
from validation import check_cfl_condition

is_stable, cfl, cfl_max = check_cfl_condition(
    dt=your_dt, dx=your_dx, diffusivity=your_D, method='explicit'
)
```

**3. Poor mesh quality warnings**

Increase minimum angle in Triangle options:
```python
B = tr.triangulate(A, 'pq33Aa')  # Try 'pq30' or higher
```

**4. Convergence issues in nonlinear iteration**

- Reduce time step
- Increase maximum iterations
- Check initial conditions

### Performance Tips

1. **Use appropriate integration rules:**
   - Linear elements (P1): `nip=3`
   - Quadratic elements (P2): `nip=7`

2. **Mesh refinement:**
   - Start coarse, refine as needed
   - Use smaller `max_area` in regions of interest

3. **Time stepping:**
   - Adaptive time stepping for nonlinear problems
   - CFL-limited for explicit methods

## Testing

Run verification tests:
```bash
python fem_shapes.py      # Test shape functions
python fem_utils.py       # Test assembly utilities
python fem_integration.py # Test quadrature rules
python validation.py      # Test validation functions
```

## References

1. Zienkiewicz, O. C., & Taylor, R. L. (2000). *The Finite Element Method, Volume 1: The Basis.* 5th ed.

2. Hughes, T. J. R. (2000). *The Finite Element Method: Linear Static and Dynamic Finite Element Analysis.*

3. Turing, A. M. (1952). The chemical basis of morphogenesis. *Philosophical Transactions of the Royal Society of London. Series B, Biological Sciences*, 237(641), 37-72.

4. Gray, P., & Scott, S. K. (1984). Autocatalytic reactions in the isothermal, continuous stirred tank reactor. *Chemical Engineering Science*, 39(6), 1087-1097.

## Migration Guide (Old → New Code)

### Shape Functions

**Old:**
```python
from shapes_tri import shapes_tri
N, dNds = shapes_tri(xi, eta)  # 2 arguments
```

**New:**
```python
from fem_shapes import shapes_tri
N, dNds = shapes_tri(xi, eta, nnodel=3)  # 3 arguments (or use default)
```

### Integration Points

**Old:**
```python
from int_points_triangle import int_points_triangle
gauss, weights = int_points_triangle(nip=3)
```

**New:**
```python
from fem_integration import ip_triangle
gauss, weights = ip_triangle(nip=3)
```

### Assembly Indices

**Old (complex nested operations):**
```python
I = np.concatenate((
    np.outer(2*EL2NOD[iel,:], np.ones(nnodel, dtype=int)).reshape(nnodel*nnodel),
    np.outer(2*EL2NOD[iel,:]+1, np.ones(nnodel, dtype=int)).reshape(nnodel*nnodel)
))
```

**New (clean and readable):**
```python
from fem_utils import create_assembly_indices_2dof
I, J = create_assembly_indices_2dof(EL2NOD, nnodel)
```

## Contributing

When adding new features:
1. Follow NumPy docstring conventions
2. Add verification/test code in `if __name__ == "__main__"` block
3. Update this README
4. Ensure backward compatibility or document breaking changes

## Contact

For questions or issues:
- Course website: [Insert URL]
- Instructor: Lars Ruepke
- GitHub Issues: [Insert URL]

---

**Last Updated:** January 2025
**Version:** 2.0 (Post-refactoring)
