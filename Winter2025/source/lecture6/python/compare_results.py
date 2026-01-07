"""
Direct comparison of Python and Numba HDF5 outputs.
"""
import numpy as np
import h5py

def compare_results(file1, file2, num_timesteps=80):
    """
    Compare temperature and heat flux results from two HDF5 files.

    The meshio XDMF format stores data as:
    - Timestep t=0: T=data2, U=data3, K=data4
    - Timestep t=1: T=data5, U=data6, K=data7
    - Pattern: T=data(2+3*t), U=data(3+3*t), K=data(4+3*t)
    """
    print(f"\n{'='*70}")
    print("VERIFICATION: Comparing Python vs Numba Results")
    print(f"{'='*70}\n")

    with h5py.File(file1, 'r') as f1, h5py.File(file2, 'r') as f2:
        max_T_diff = 0.0
        max_U_diff = 0.0
        max_K_diff = 0.0

        mean_T_diff = 0.0
        mean_U_diff = 0.0

        for t in range(num_timesteps):
            # Calculate data indices for this timestep
            T_idx = 2 + 3 * t
            U_idx = 3 + 3 * t
            K_idx = 4 + 3 * t

            # Read temperature field
            T1 = f1[f'data{T_idx}'][:]
            T2 = f2[f'data{T_idx}'][:]
            T_diff = np.abs(T1 - T2)
            max_T_diff = max(max_T_diff, np.max(T_diff))
            mean_T_diff += np.mean(T_diff)

            # Read heat flux
            U1 = f1[f'data{U_idx}'][:]
            U2 = f2[f'data{U_idx}'][:]
            U_diff = np.abs(U1 - U2)
            max_U_diff = max(max_U_diff, np.max(U_diff))
            mean_U_diff += np.mean(U_diff)

            # Read conductivity (should be identical - not computed)
            K1 = f1[f'data{K_idx}'][:]
            K2 = f2[f'data{K_idx}'][:]
            K_diff = np.abs(K1 - K2)
            max_K_diff = max(max_K_diff, np.max(K_diff))

            # Print progress every 20 timesteps
            if (t + 1) % 20 == 0 or t == 0 or t == num_timesteps - 1:
                print(f"  Timestep {t:3d}: "
                      f"max ΔT = {np.max(T_diff):.2e}, "
                      f"max ΔU = {np.max(U_diff):.2e}, "
                      f"mean ΔT = {np.mean(T_diff):.2e}")

        mean_T_diff /= num_timesteps
        mean_U_diff /= num_timesteps

    print(f"\n{'='*70}")
    print("VERIFICATION SUMMARY")
    print(f"{'='*70}")
    print(f"Temperature field:")
    print(f"  Maximum absolute difference: {max_T_diff:.3e}")
    print(f"  Mean absolute difference:    {mean_T_diff:.3e}")
    print(f"\nHeat flux field:")
    print(f"  Maximum absolute difference: {max_U_diff:.3e}")
    print(f"  Mean absolute difference:    {mean_U_diff:.3e}")
    print(f"\nConductivity field:")
    print(f"  Maximum absolute difference: {max_K_diff:.3e} (should be 0)")

    # Determine pass/fail
    tolerance = 1e-10
    passed = (max_T_diff < tolerance and max_U_diff < tolerance and max_K_diff < tolerance)

    print(f"\n{'─'*70}")
    if passed:
        print("✓ VERIFICATION PASSED")
        print("  Results are numerically identical within tolerance (< 1e-10)")
    else:
        print("✗ VERIFICATION FAILED")
        print(f"  Differences exceed tolerance ({tolerance:.2e})")

    print(f"{'='*70}\n")

    return passed

if __name__ == "__main__":
    success = compare_results('transient_python.h5', 'transient_numba.h5', num_timesteps=80)
    exit(0 if success else 1)
