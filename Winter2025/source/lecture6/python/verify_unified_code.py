#!/usr/bin/env python3
"""
Verification script for unified Python/Numba code.

Tests that the refactored v3 code (with identical Python/Numba implementations)
produces the same results in both modes.
"""

import numpy as np
import h5py
import subprocess
import sys
from pathlib import Path

def modify_use_numba(filename: str, use_numba: bool) -> None:
    """Modify the USE_NUMBA flag in the source file."""
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Find and modify the USE_NUMBA line
    for i, line in enumerate(lines):
        if line.strip().startswith('USE_NUMBA ='):
            if use_numba:
                lines[i] = 'USE_NUMBA = True\n'
            else:
                lines[i] = 'USE_NUMBA = False\n'
            break

    with open(filename, 'w') as f:
        f.writelines(lines)

def run_simulation(mode: str, output_file: str) -> bool:
    """Run simulation in specified mode."""
    print(f"\n{'='*60}")
    print(f"Running simulation in {mode} mode...")
    print(f"{'='*60}")

    # Modify the code to use the specified mode
    use_numba = (mode == "NUMBA")
    modify_use_numba('2d_fem_transient_triangle_v3.py', use_numba)

    # Run the simulation
    result = subprocess.run(
        [sys.executable, '2d_fem_transient_triangle_v3.py'],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print(f"ERROR: Simulation failed in {mode} mode")
        print(result.stderr)
        return False

    # Rename output file
    if Path('transient.h5').exists():
        Path('transient.h5').rename(output_file)
        print(f"✓ Output saved to {output_file}")
    else:
        print(f"ERROR: Output file not created")
        return False

    return True

def compare_results(file1: str, file2: str) -> bool:
    """Compare results from two HDF5 files."""
    print(f"\n{'='*60}")
    print("Comparing results...")
    print(f"{'='*60}")

    with h5py.File(file1, 'r') as f1, h5py.File(file2, 'r') as f2:
        # Get all data keys (meshio format: data0, data1, data2, ...)
        keys1 = sorted([k for k in f1.keys() if k.startswith('data')],
                       key=lambda x: int(x[4:]))
        keys2 = sorted([k for k in f2.keys() if k.startswith('data')],
                       key=lambda x: int(x[4:]))

        # Check same number of datasets
        if len(keys1) != len(keys2):
            print(f"ERROR: Different number of datasets!")
            print(f"  Python: {len(keys1)}")
            print(f"  Numba:  {len(keys2)}")
            return False

        print(f"Number of datasets: {len(keys1)}")

        # In meshio format: data0, data1, data2 = x, y, T for first timestep
        # Then data3, data4, data5 = Qx, Qy, T for same timestep
        # Pattern repeats with output_freq

        max_diff_T = 0.0
        max_diff_Qx = 0.0
        max_diff_Qy = 0.0

        # Compare each dataset
        n_compared = 0
        for key1, key2 in zip(keys1, keys2):
            data1 = f1[key1][:]
            data2 = f2[key2][:]

            # Check shapes match
            if data1.shape != data2.shape:
                print(f"ERROR: Shape mismatch in {key1}")
                print(f"  Python: {data1.shape}")
                print(f"  Numba:  {data2.shape}")
                return False

            # Compute difference
            diff = np.abs(data1 - data2)
            max_diff = np.max(diff)

            # Categorize by variable name (stored in attributes)
            var_name = f1[key1].attrs.get('name', 'unknown').decode() if isinstance(f1[key1].attrs.get('name', 'unknown'), bytes) else f1[key1].attrs.get('name', 'unknown')

            if var_name == 'T':
                max_diff_T = max(max_diff_T, max_diff)
            elif var_name == 'Q_x':
                max_diff_Qx = max(max_diff_Qx, max_diff)
            elif var_name == 'Q_y':
                max_diff_Qy = max(max_diff_Qy, max_diff)

            n_compared += 1

            # Warn if difference is large
            tol = 1e-12
            if max_diff > tol and var_name not in ['x', 'y']:
                print(f"\nWARNING: {key1} ({var_name}) has differences > {tol}")
                print(f"  Max diff = {max_diff:.3e}")

        # Summary
        print(f"\n{'='*60}")
        print("RESULTS SUMMARY")
        print(f"{'='*60}")
        print(f"Compared {n_compared} datasets")
        print(f"Maximum differences:")
        print(f"  Temperature: {max_diff_T:.3e}")
        print(f"  Heat flux X: {max_diff_Qx:.3e}")
        print(f"  Heat flux Y: {max_diff_Qy:.3e}")

        if max_diff_T < 1e-12 and max_diff_Qx < 1e-12 and max_diff_Qy < 1e-12:
            print(f"\n✓ VERIFICATION PASSED")
            print(f"  Results are identical within machine precision (< 1e-12)")
            print(f"  Python and Numba modes produce the same output!")
            return True
        else:
            print(f"\n✗ VERIFICATION FAILED")
            print(f"  Differences exceed machine precision")
            return False

def main():
    """Main verification routine."""
    print("="*60)
    print("UNIFIED CODE VERIFICATION")
    print("="*60)
    print("Testing that unified Python/Numba implementation")
    print("produces identical results in both modes.")
    print("="*60)

    # Output files
    python_output = "transient_python.h5"
    numba_output = "transient_numba.h5"

    # Run Python mode
    if not run_simulation("PYTHON", python_output):
        print("\n✗ Python mode failed")
        return False

    # Run Numba mode
    if not run_simulation("NUMBA", numba_output):
        print("\n✗ Numba mode failed")
        return False

    # Compare results
    results_match = compare_results(python_output, numba_output)

    # Restore default (Numba on)
    modify_use_numba('2d_fem_transient_triangle_v3.py', True)

    return results_match

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
