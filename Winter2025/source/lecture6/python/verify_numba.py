"""
Verification script to compare Python and Numba implementations.
Runs both modes and compares temperature and heat flux results.
"""
import numpy as np
import h5py
import subprocess
import sys

def modify_use_numba(filename, use_numba_value):
    """Modify the USE_NUMBA flag in the source file."""
    with open(filename, 'r') as f:
        content = f.read()

    # Replace the USE_NUMBA line
    if use_numba_value:
        content = content.replace('USE_NUMBA = False', 'USE_NUMBA = True')
    else:
        content = content.replace('USE_NUMBA = True', 'USE_NUMBA = False')

    with open(filename, 'w') as f:
        f.write(content)

def modify_output_filename(filename, output_name):
    """Modify the output filename in the source file."""
    with open(filename, 'r') as f:
        content = f.read()

    # Replace the output filename
    content = content.replace("'transient.xmf'", f"'{output_name}.xmf'")

    with open(filename, 'w') as f:
        f.write(content)

def run_simulation(mode_name):
    """Run simulation and return status."""
    print(f"\n{'='*70}")
    print(f"Running simulation in {mode_name} mode...")
    print(f"{'='*70}")

    result = subprocess.run(
        ['python', '2d_fem_transient_triangle_v3.py'],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print(f"Error running {mode_name} mode:")
        print(result.stderr)
        return False

    # Print last few lines of output
    lines = result.stdout.strip().split('\n')
    for line in lines[-15:]:
        print(line)

    return True

def compare_hdf5_files(file1, file2, tolerance=1e-10):
    """Compare two HDF5 files and report differences."""
    print(f"\n{'='*70}")
    print("COMPARING RESULTS")
    print(f"{'='*70}")

    with h5py.File(file1, 'r') as f1, h5py.File(file2, 'r') as f2:
        # Get list of time steps (meshio uses 'data' prefix)
        all_keys = list(f1.keys())
        timesteps = sorted([key for key in all_keys if key.startswith('data')],
                          key=lambda x: int(x.replace('data', '')))

        print(f"\nFound {len(timesteps)} datasets to compare")

        max_T_diff = 0.0
        max_U_diff = 0.0
        max_K_diff = 0.0

        for timestep in timesteps:
            step_num = int(timestep.split('_')[1])

            # Compare temperature field (T)
            T1 = f1[timestep]['T'][:]
            T2 = f2[timestep]['T'][:]
            T_diff = np.abs(T1 - T2)
            max_T_diff = max(max_T_diff, np.max(T_diff))

            # Compare heat flux (U)
            U1 = f1[timestep]['U'][:]
            U2 = f2[timestep]['U'][:]
            U_diff = np.abs(U1 - U2)
            max_U_diff = max(max_U_diff, np.max(U_diff))

            # Compare conductivity (K)
            K1 = f1[timestep]['K'][:]
            K2 = f2[timestep]['K'][:]
            K_diff = np.abs(K1 - K2)
            max_K_diff = max(max_K_diff, np.max(K_diff))

            # Report progress every 20 steps
            if step_num % 20 == 0 or step_num == len(timesteps):
                print(f"  Step {step_num:3d}: max T diff = {np.max(T_diff):.2e}, "
                      f"max U diff = {np.max(U_diff):.2e}, "
                      f"max K diff = {np.max(K_diff):.2e}")

        print(f"\n{'='*70}")
        print("VERIFICATION SUMMARY")
        print(f"{'='*70}")
        print(f"Maximum temperature difference:  {max_T_diff:.2e}")
        print(f"Maximum heat flux difference:    {max_U_diff:.2e}")
        print(f"Maximum conductivity difference: {max_K_diff:.2e}")
        print(f"Tolerance threshold:             {tolerance:.2e}")

        # Check if differences are within tolerance
        if max_T_diff < tolerance and max_U_diff < tolerance and max_K_diff < tolerance:
            print(f"\n✓ VERIFICATION PASSED: Results are identical within tolerance!")
            print(f"{'='*70}\n")
            return True
        else:
            print(f"\n✗ VERIFICATION FAILED: Differences exceed tolerance!")
            print(f"{'='*70}\n")
            return False

def main():
    source_file = '2d_fem_transient_triangle_v3.py'

    # Read original content to restore later
    with open(source_file, 'r') as f:
        original_content = f.read()

    try:
        # Run Python mode
        print("\n" + "="*70)
        print("STEP 1: PYTHON MODE")
        print("="*70)
        modify_use_numba(source_file, False)
        modify_output_filename(source_file, 'transient_python')
        if not run_simulation("PYTHON"):
            return 1

        # Restore original and run Numba mode
        with open(source_file, 'w') as f:
            f.write(original_content)

        print("\n" + "="*70)
        print("STEP 2: NUMBA MODE")
        print("="*70)
        modify_use_numba(source_file, True)
        modify_output_filename(source_file, 'transient_numba')
        if not run_simulation("NUMBA"):
            return 1

        # Compare results
        success = compare_hdf5_files('transient_python.h5', 'transient_numba.h5')

        return 0 if success else 1

    finally:
        # Restore original file
        with open(source_file, 'w') as f:
            f.write(original_content)
        print("Original file restored.")

if __name__ == "__main__":
    sys.exit(main())
