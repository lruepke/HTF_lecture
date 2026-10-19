# Numba vs Python Assembly Loops: Key Differences

## Overview

The v3 code contains **two parallel implementations** of the element assembly loop:
- **Python version**: Original implementation, interpreted
- **Numba version**: Compiled to machine code via JIT

Both produce **identical results** but with very different performance characteristics.

## Side-by-Side Comparison

### Python Version (Inside `solve_2d_temperature_fem`)

```python
# ========== PYTHON PATH ==========
# Storage
Rhs_all = np.zeros(nnod)
I = np.zeros((nel, nnodel*nnodel))
J = np.zeros((nel, nnodel*nnodel))
K = np.zeros((nel, nnodel*nnodel))

# Element assembly loop
for iel in range(0, nel):
    # Access mesh data through object attributes
    ECOORD = np.take(mesh.GCOORD, mesh.EL2NOD[iel, :], axis=0)
    Ael = np.zeros((nnodel, nnodel))
    Rhs_el = np.zeros(nnodel)

    for ip in range(0, integration.nip):
        # Get pre-computed shape functions from object
        N = integration.N_all[ip, :]
        dNds = integration.dNds_all[ip, :, :]

        # Jacobian calculations
        Jac = np.matmul(dNds, ECOORD)
        invJ = np.linalg.inv(Jac)
        detJ = np.linalg.det(Jac)
        dNdx = np.matmul(invJ, dNds)

        # Element matrices
        Me = np.outer(N, N)
        Ke = dNdx.T @ dNdx

        # Assemble element matrix (accessing object attributes)
        Ael += (material.rho*material.cp*Me + time_params.dt*Kel[iel]*Ke) * detJ * integration.weights[ip]

        # RHS (accessing mesh connectivity through object)
        Rhs_el += material.rho*material.cp * (Me @ T[mesh.EL2NOD[iel, :]]) * detJ * integration.weights[ip]

    # Assemble global arrays using numpy broadcasting
    I[iel, :] = (mesh.EL2NOD[iel, :]*np.ones((nnodel, 1), dtype=int)).T.reshape(nnodel*nnodel)
    J[iel, :] = (mesh.EL2NOD[iel, :]*np.ones((nnodel, 1), dtype=int)).reshape(nnodel*nnodel)
    K[iel, :] = Ael.reshape(nnodel*nnodel)

    # Accumulate RHS (accessing mesh connectivity through object)
    Rhs_all[mesh.EL2NOD[iel, :]] += Rhs_el
```

### Numba Version (Separate Function)

```python
@njit(cache=True)
def assemble_system_numba(
    GCOORD: np.ndarray,      # All data passed as parameters
    EL2NOD: np.ndarray,
    Kel: np.ndarray,
    T: np.ndarray,
    N_all: np.ndarray,       # Pre-computed shape functions
    dNds_all: np.ndarray,
    weights: np.ndarray,
    rho: float,              # Scalars passed directly
    cp: float,
    dt: float
) -> tuple:
    """Numba-optimized element assembly."""
    nel = EL2NOD.shape[0]
    nnodel = EL2NOD.shape[1]
    nnod = GCOORD.shape[0]
    nip = N_all.shape[0]

    # Storage
    Rhs_all = np.zeros(nnod)
    I = np.zeros((nel, nnodel * nnodel))
    J = np.zeros((nel, nnodel * nnodel))
    K = np.zeros((nel, nnodel * nnodel))

    for iel in range(nel):
        # Direct array indexing (no np.take)
        ECOORD = GCOORD[EL2NOD[iel, :], :]
        Ael = np.zeros((nnodel, nnodel))
        Rhs_el = np.zeros(nnodel)

        for ip in range(nip):
            # Direct array access
            N = N_all[ip, :]
            dNds = dNds_all[ip, :, :]

            # Jacobian calculations (same operations)
            Jac = dNds @ ECOORD
            invJ = np.linalg.inv(Jac)
            detJ = np.linalg.det(Jac)
            dNdx = invJ @ dNds

            # Element matrices (same operations)
            Me = np.outer(N, N)
            Ke_local = dNdx.T @ dNdx

            # Assemble (scalars passed as parameters, not attributes)
            Ael += (rho * cp * Me + dt * Kel[iel] * Ke_local) * detJ * weights[ip]

            # RHS (direct indexing)
            T_el = T[EL2NOD[iel, :]]
            Rhs_el += rho * cp * (Me @ T_el) * detJ * weights[ip]

        # Store coefficients - explicit loops instead of broadcasting
        for i in range(nnodel):
            for j in range(nnodel):
                idx = i * nnodel + j
                I[iel, idx] = EL2NOD[iel, i]
                J[iel, idx] = EL2NOD[iel, j]
                K[iel, idx] = Ael[i, j]

        # Add to global RHS - explicit loop
        for i in range(nnodel):
            Rhs_all[EL2NOD[iel, i]] += Rhs_el[i]

    return I, J, K, Rhs_all
```

## Key Differences

### 1. Function Scope and Data Access

**Python:**
- Embedded in main solver function
- Accesses data through objects: `mesh.GCOORD`, `material.rho`, `integration.weights`
- Can use any Python feature

**Numba:**
- Separate standalone function
- All data passed as **explicit parameters**
- No object attribute access (numba can't optimize that)
- Limited to numba-supported Python subset

### 2. Array Operations

**Python:**
```python
# Uses numpy convenience functions
ECOORD = np.take(mesh.GCOORD, mesh.EL2NOD[iel, :], axis=0)

# Fancy indexing with broadcasting
I[iel, :] = (mesh.EL2NOD[iel, :]*np.ones((nnodel, 1), dtype=int)).T.reshape(nnodel*nnodel)
```

**Numba:**
```python
# Direct indexing (faster in compiled code)
ECOORD = GCOORD[EL2NOD[iel, :], :]

# Explicit loops (compiler can optimize these)
for i in range(nnodel):
    for j in range(nnodel):
        idx = i * nnodel + j
        I[iel, idx] = EL2NOD[iel, i]
```

**Why the difference?**
- In interpreted Python, vectorized operations are faster (implemented in C)
- In compiled code, explicit loops are often faster (compiler can optimize)
- Numba generates SIMD instructions for simple loops

### 3. Function Call Pattern

**Python:**
```python
# Called directly in solve function
if USE_NUMBA and NUMBA_AVAILABLE:
    # ... numba path ...
else:
    # Python loop runs here inline
    for iel in range(0, nel):
        # ...
```

**Numba:**
```python
# Called as a separate function
I, J, K, Rhs_all = assemble_system_numba(
    mesh.GCOORD,           # Extract arrays from objects
    mesh.EL2NOD,
    Kel,
    T,
    integration.N_all,     # Pass pre-computed data
    integration.dNds_all,
    integration.weights,
    material.rho,          # Pass scalars
    material.cp,
    time_params.dt
)
```

### 4. Return Values

**Python:**
- Modifies arrays in-place (`I`, `J`, `K`, `Rhs_all` already exist)
- No explicit return of assembly results

**Numba:**
- Returns assembled arrays as tuple
- Cleaner functional style
- Caller receives: `I, J, K, Rhs_all = assemble_system_numba(...)`

### 5. Type Information

**Python:**
```python
# No type constraints
for iel in range(0, nel):  # nel inferred at runtime
    ECOORD = np.take(...)  # shape unknown until runtime
```

**Numba:**
```python
@njit(cache=True)  # Types inferred on first call, then cached
def assemble_system_numba(
    GCOORD: np.ndarray,  # Type hints for documentation only
    ...                  # Numba infers actual types from usage
) -> tuple:
```

Numba **infers types** on first call:
- First call: slow (compiles function based on input types)
- Subsequent calls: fast (uses cached compiled version)

## Performance Impact

### Why Numba is Faster

1. **Compilation to Machine Code**
   - Python: Interpreted line-by-line
   - Numba: Compiled to native CPU instructions

2. **Loop Optimization**
   - Python: Loop overhead every iteration
   - Numba: Loops unrolled, pipelined by compiler

3. **Type Specialization**
   - Python: Type checks at runtime
   - Numba: Types known at compile time, no checks needed

4. **SIMD Vectorization**
   - Python: Limited vectorization (only in numpy functions)
   - Numba: Compiler generates SIMD instructions for loops

5. **Memory Access**
   - Python: Pointer chasing through object attributes
   - Numba: Direct memory access, better cache utilization

### The Cost of Numba

**First call is slow** (JIT compilation):
```
Time step 1 (with compilation): 3078 ms  ← Compiling!
Time step 2 (cached):              8.5 ms  ← Fast
```

This is why we use `cache=True`:
- Compiled version saved to disk
- Subsequent program runs skip compilation
- Only recompiles if code changes

## Trade-offs

| Aspect | Python | Numba |
|--------|--------|-------|
| **Development speed** | ✓ Fast (any Python) | Slower (subset only) |
| **Debugging** | ✓ Easy | Harder (compiled) |
| **Flexibility** | ✓ Full Python | Limited features |
| **Runtime speed** | Slow (interpreted) | ✓ Fast (compiled) |
| **First run time** | ✓ Fast | Slow (compilation) |
| **Dependencies** | ✓ None | Requires numba |

## Why Keep Both?

The v3 code maintains **both implementations** for good reasons:

1. **Educational**: Students see both approaches
2. **Debugging**: Python version easier to debug
3. **Validation**: Compare results to ensure correctness
4. **Portability**: Works without numba installation
5. **Flexibility**: Easy to toggle for development vs production

## When to Use Each

**Use Python version when:**
- Developing/debugging new features
- Making small changes and testing
- Teaching students about algorithms
- Numba not available

**Use Numba version when:**
- Running production simulations
- Need maximum performance
- Code is stable and tested
- Have numba installed

## Example: Switching Between Modes

Simply change one line at the top of the file:

```python
# ============================================================================
# PERFORMANCE CONFIGURATION
# ============================================================================
USE_NUMBA = True   # Fast: 14x speedup
# USE_NUMBA = False  # Slow: easier to debug
```

Everything else remains the same!

## Summary

The core **algorithm is identical** in both versions:
- Same mathematical operations
- Same element assembly logic
- Same integration scheme
- Same sparse matrix format

The **implementation differs** to meet numba's requirements:
- Separate function (not inline)
- Explicit parameters (not object attributes)
- Simpler indexing (not fancy numpy)
- Explicit loops (not broadcasting)

**Result**: 14x speedup while keeping both implementations available for teaching, debugging, and validation.
