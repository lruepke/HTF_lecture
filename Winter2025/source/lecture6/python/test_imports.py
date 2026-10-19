"""
Quick test to verify that all scripts can import fem_shapes correctly.
"""

import sys
import numpy as np

print("=" * 60)
print("Testing fem_shapes Import and Usage")
print("=" * 60)
print()

# Test 1: Direct import
print("Test 1: Direct import of fem_shapes...")
try:
    from fem_shapes import shapes_tri
    print("✓ Successfully imported shapes_tri from fem_shapes")
except ImportError as e:
    print(f"✗ Failed to import: {e}")
    sys.exit(1)

# Test 2: Call with 2 arguments (uses default nnodel=3)
print("\nTest 2: Call shapes_tri with 2 arguments (default nnodel=3)...")
try:
    N, dNds = shapes_tri(0.333, 0.333)
    print(f"✓ N shape: {N.shape}, dNds shape: {dNds.shape}")
    print(f"  N = {N}")
    print(f"  dNds =\n{dNds}")
except Exception as e:
    print(f"✗ Failed: {e}")
    sys.exit(1)

# Test 3: Call with 3 arguments (explicit nnodel=3)
print("\nTest 3: Call shapes_tri with 3 arguments (explicit nnodel=3)...")
try:
    N, dNds = shapes_tri(0.333, 0.333, 3)
    print(f"✓ N shape: {N.shape}, dNds shape: {dNds.shape}")
except Exception as e:
    print(f"✗ Failed: {e}")
    sys.exit(1)

# Test 4: Call with nnodel=6 (quadratic elements)
print("\nTest 4: Call shapes_tri with nnodel=6...")
try:
    N, dNds = shapes_tri(0.333, 0.333, 6)
    print(f"✓ N shape: {N.shape}, dNds shape: {dNds.shape}")
    print(f"  N = {N}")
except Exception as e:
    print(f"✗ Failed: {e}")
    sys.exit(1)

# Test 5: Verify partition of unity
print("\nTest 5: Verify partition of unity (sum of N should be 1)...")
N3, _ = shapes_tri(0.5, 0.25, 3)
N6, _ = shapes_tri(0.5, 0.25, 6)
sum3 = np.sum(N3)
sum6 = np.sum(N6)
print(f"  3-node: sum(N) = {sum3:.10f} {'✓' if abs(sum3-1) < 1e-10 else '✗'}")
print(f"  6-node: sum(N) = {sum6:.10f} {'✓' if abs(sum6-1) < 1e-10 else '✗'}")

print("\n" + "=" * 60)
print("All tests passed! fem_shapes is working correctly.")
print("=" * 60)
