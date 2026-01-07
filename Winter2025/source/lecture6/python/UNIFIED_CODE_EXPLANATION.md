# Unified Python/Numba Implementation

## The Key Insight: IDENTICAL Code

The v3 code now uses a **unified approach** where the Python and Numba implementations are **literally the same code**. The only difference is whether that code gets JIT-compiled by numba or runs as regular Python.

## How It Works

### 1. Write Functions Once (Numba-Compatible Style)

Both assembly functions are written **once**, in a style that works with or without numba:

```python
def assemble_system_optimized(
    GCOORD: np.ndarray,
    EL2NOD: np.ndarray,
    Kel: np.ndarray,
    T: np.ndarray,
    N_all: np.ndarray,
    dNds_all: np.ndarray,
    weights: np.ndarray,
    rho: float,
    cp: float,
    dt: float
) -> tuple:
    """
    Element assembly loop.

    This function works as pure Python OR can be JIT-compiled by numba.
    The code is IDENTICAL in both cases - only execution speed differs.
    """
    nel = EL2NOD.shape[0]
    nnodel = EL2NOD.shape[1]
    # ... rest of implementation ...

    for iel in range(nel):
        ECOORD = GCOORD[EL2NOD[iel, :], :]
        # ... element assembly ...

    return I, J, K, Rhs_all
```

Key characteristics of numba-compatible code:
- **No object attribute access** (`mesh.GCOORD` → `GCOORD` as parameter)
- **Simple indexing** (`array[indices]` not `np.take()`)
- **Explicit loops** (instead of fancy broadcasting)
- **Standard numpy functions** that numba supports

### 2. Conditionally Apply JIT Compilation

After defining the functions, we **optionally** apply numba's JIT compiler:

```python
# ============================================================================
# APPLY NUMBA JIT COMPILATION (if enabled)
# ============================================================================
# The functions above are written in a numba-compatible style.
# If USE_NUMBA is True, we compile them to machine code for 10-20x speedup.
# If False, they run as normal Python functions.
# THE CODE IS IDENTICAL - only the execution speed differs!

if USE_NUMBA and NUMBA_AVAILABLE:
    print("Applying numba JIT compilation to assembly functions...")
    assemble_system_optimized = njit(cache=True)(assemble_system_optimized)
    compute_heat_flux_optimized = njit(cache=True)(compute_heat_flux_optimized)
```

**What this does:**
- If `USE_NUMBA = True`: Functions get wrapped with `@njit` → compiled to machine code
- If `USE_NUMBA = False`: Functions remain as regular Python → no compilation

### 3. Call the Same Function Regardless

The solver just calls the function - it doesn't need to know whether it's compiled or not:

```python
def solve_2d_temperature_fem(...):
    # ...

    # Element assembly
    # Same code runs whether numba is on or off - only execution speed differs!
    t_start_assembly = time.perf_counter()
    I, J, K, Rhs_all = assemble_system_optimized(
        mesh.GCOORD,
        mesh.EL2NOD,
        Kel,
        T,
        integration.N_all,
        integration.dNds_all,
        integration.weights,
        material.rho,
        material.cp,
        time_params.dt
    )
    t_assembly = time.perf_counter() - t_start_assembly

    # ... rest of solver ...
```

**No conditional logic needed!** The same function call works either way.

## Visual Comparison

### OLD Approach (Dual Code Paths)

```python
# Two completely separate implementations:

if USE_NUMBA and NUMBA_AVAILABLE:
    # ========== NUMBA PATH ==========
    @njit(cache=True)
    def assemble_system_numba(...):
        for iel in range(nel):
            ECOORD = GCOORD[EL2NOD[iel, :], :]
            # ... assembly code ...

    # Call numba version
    I, J, K, Rhs_all = assemble_system_numba(...)

else:
    # ========== PYTHON PATH ==========
    for iel in range(nel):
        ECOORD = np.take(mesh.GCOORD, mesh.EL2NOD[iel, :], axis=0)
        # ... assembly code (slightly different) ...
```

**Problems:**
- Code duplication
- Two versions to maintain
- Risk of divergence
- Harder to teach

### NEW Approach (Unified Code)

```python
# Single implementation, conditionally compiled:

def assemble_system_optimized(...):
    """Works as Python OR gets compiled by numba"""
    for iel in range(nel):
        ECOORD = GCOORD[EL2NOD[iel, :], :]
        # ... assembly code ...
    return I, J, K, Rhs_all

# Optionally compile
if USE_NUMBA and NUMBA_AVAILABLE:
    assemble_system_optimized = njit(cache=True)(assemble_system_optimized)

# Call the same function regardless
I, J, K, Rhs_all = assemble_system_optimized(...)
```

**Benefits:**
- **Zero code duplication** - write once, run both ways
- **Same algorithm** - impossible to diverge
- **Easy to maintain** - only one version to update
- **Better for teaching** - shows that numba is just compiling Python!

## Performance Comparison

### Pure Python Mode (`USE_NUMBA = False`)

```
Time step 10/80 - Performance breakdown (Mode: PYTHON):
  Element assembly:      106.58 ms (97.5%)
  Linear solve:            2.42 ms (2.2%)
  ──────────────────────────────────────
  Total time:            109.29 ms
```

### Numba Mode (`USE_NUMBA = True`)

```
Time step 10/80 - Performance breakdown (Mode: NUMBA):
  Element assembly:        7.57 ms (77.0%)
  Linear solve:            2.07 ms (21.1%)
  ──────────────────────────────────────
  Total time:              9.83 ms
```

**Same code, 11x speedup!**

## Verification Results

The unified code has been verified to produce **bit-for-bit identical results** in both Python and Numba modes:

```
============================================================
RESULTS SUMMARY
============================================================
Compared 29 datasets
Maximum differences:
  Temperature: 0.000e+00
  Heat flux X: 0.000e+00
  Heat flux Y: 0.000e+00

✓ VERIFICATION PASSED
  Results are identical within machine precision (< 1e-12)
  Python and Numba modes produce the same output!
```

This confirms that:
- The refactoring was successful
- Both modes use the same algorithm
- No numerical differences introduced
- Perfect agreement across all timesteps and variables

See `UNIFIED_CODE_VERIFICATION.md` for detailed test results.

## What's Different From Before?

### Before (v3 original):
```python
# Separate numba function
if USE_NUMBA and NUMBA_AVAILABLE:
    @njit(cache=True)
    def assemble_system_numba(...):
        # implementation

# Separate Python loop in solver
else:
    for iel in range(nel):
        # slightly different implementation
```

### After (v3 unified):
```python
# One function definition
def assemble_system_optimized(...):
    # implementation

# Conditional compilation
if USE_NUMBA and NUMBA_AVAILABLE:
    assemble_system_optimized = njit(cache=True)(assemble_system_optimized)

# Always call the same function
I, J, K, Rhs_all = assemble_system_optimized(...)
```

## Educational Value

This unified approach teaches several important concepts:

1. **Numba compiles Python** - it's not a different language!
2. **Performance through compilation** - same algorithm, different execution
3. **Code reuse** - write once, run multiple ways
4. **Decorators** - how `@njit` modifies function behavior
5. **Abstraction** - the caller doesn't need to know implementation details

## How to Toggle

Simply change one line at the top of the file:

```python
# Fast: 11x speedup, requires numba
USE_NUMBA = True

# Slow: easier to debug, no dependencies
# USE_NUMBA = False
```

Everything else is automatic!

## Code Style Guidelines for Numba Compatibility

To write Python code that can be compiled by numba:

### ✅ DO:
- Pass arrays and scalars as function parameters
- Use simple array indexing: `array[i, j]`
- Use explicit loops where needed
- Use standard numpy functions (inv, det, outer, @, etc.)
- Keep functions pure (no global state)

### ❌ DON'T:
- Access object attributes inside loops: `mesh.GCOORD`
- Use fancy numpy helpers: `np.take()`, complex broadcasting
- Use Python features numba doesn't support (dicts, classes, etc.)
- Modify global variables
- Use non-numpy libraries

## Summary

**Old approach**: Two separate code paths, duplicated logic

**New approach**: One code path, conditionally compiled

**Result**:
- Simpler code
- Easier to maintain
- Impossible for implementations to diverge
- Same 11x performance boost
- Better teaching tool - shows numba is just compiling Python!

The key insight: **Write Python code in a numba-friendly style, then optionally compile it.** No separate implementations needed!
