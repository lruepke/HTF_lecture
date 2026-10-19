"""
Mesh Generation Utilities for Triangle Meshes

This module provides helper functions for creating triangle meshes using
the Triangle library (https://www.cs.cmu.edu/~quake/triangle.html).

These utilities define geometric domains (boxes, circles) that can be
combined to create complex mesh geometries for finite element analysis.

Author: HTF Lecture Series
Date: 2025
"""

import numpy as np


def make_box(vertices, segments, regions, x, y, w, h, attribute, max_area=0.005):
    """
    Add a rectangular domain to mesh definition.

    This function adds a rectangular box to the mesh specification by
    appending vertices, segments, and a region marker. The box is defined
    by its bottom-left corner (x, y) and dimensions (w, h).

    Parameters
    ----------
    vertices : list
        List of [x, y] vertex coordinates.
        Modified in-place by adding 4 corner vertices.
    segments : list
        List of (i, j) segment endpoint indices.
        Modified in-place by adding 4 edge segments.
    regions : list
        List of region specifications [x, y, attribute, max_area].
        Modified in-place by adding 1 region marker.
    x, y : float
        Coordinates of bottom-left corner [m]
    w, h : float
        Width and height of rectangle [m]
    attribute : int
        Region attribute/material ID (positive integer).
        Used to assign different material properties to different regions.
    max_area : float, optional
        Maximum triangle area in this region [m²].
        Smaller values create finer meshes. Default is 0.005.

    Returns
    -------
    None
        Modifies vertices, segments, and regions lists in-place.

    Notes
    -----
    The box is created with 4 vertices in counter-clockwise order:
        Vertex 0: (x, y)       - bottom-left
        Vertex 1: (x+w, y)     - bottom-right
        Vertex 2: (x+w, y+h)   - top-right
        Vertex 3: (x, y+h)     - top-left

    The region marker is placed slightly inside the box at (x+0.01*w, y+0.01*h)
    to ensure it's inside the domain for the Triangle library.

    Typically used to define the outer domain boundary, but can also be used
    for internal rectangular features.

    Examples
    --------
    Create a 1m x 1m square domain:

    >>> vertices, segments, regions = [], [], []
    >>> make_box(vertices, segments, regions, 0, 0, 1, 1, attribute=1)
    >>> print(f"Created {len(vertices)} vertices, {len(segments)} segments")
    Created 4 vertices, 4 segments

    Create domain with finer mesh:

    >>> vertices, segments, regions = [], [], []
    >>> make_box(vertices, segments, regions, 0, 0, 1, 1,
    ...          attribute=1, max_area=0.001)
    >>> # Will create more triangles due to smaller max_area

    See Also
    --------
    make_inclusion : Add circular inclusion to mesh
    create_triangle_mesh : Complete mesh creation workflow
    """
    # Get current number of vertices (for indexing new vertices)
    i = len(vertices)

    # Add 4 corner vertices in counter-clockwise order
    vertices.extend([
        [x,     y],      # Bottom-left
        [x + w, y],      # Bottom-right
        [x + w, y + h],  # Top-right
        [x,     y + h]   # Top-left
    ])

    # Add 4 edge segments (connecting consecutive vertices)
    segments.extend([
        (i + 0, i + 1),  # Bottom edge
        (i + 1, i + 2),  # Right edge
        (i + 2, i + 3),  # Top edge
        (i + 3, i + 0)   # Left edge
    ])

    # Add region marker (inside the box, with attribute and max area)
    # Region format: [x_inside, y_inside, attribute, max_area]
    # Place marker slightly inside to avoid boundary issues
    regions.append([x + 0.01*w, y + 0.01*h, attribute, max_area])


def make_inclusion(vertices, segments, regions, center_x, center_y, radius,
                   n_points, attribute, max_area=0.001):
    """
    Add a circular inclusion to mesh definition.

    This function adds a circular region to the mesh specification, commonly
    used to represent circular heterogeneities like low-conductivity inclusions
    in heat transfer problems or different material phases.

    Parameters
    ----------
    vertices : list
        List of [x, y] vertex coordinates.
        Modified in-place by adding n_points vertices around circle perimeter.
    segments : list
        List of (i, j) segment endpoint indices.
        Modified in-place by adding n_points segments forming closed loop.
    regions : list
        List of region specifications [x, y, attribute, max_area].
        Modified in-place by adding 1 region marker at circle center.
    center_x, center_y : float
        Circle center coordinates [m]
    radius : float
        Circle radius [m]
    n_points : int
        Number of points on circle perimeter.
        Higher values create smoother circular boundary (typical: 20-40).
    attribute : int
        Region attribute/material ID (positive integer).
        Typically different from surrounding domain to distinguish materials.
    max_area : float, optional
        Maximum triangle area inside this region [m²].
        Smaller values create finer mesh. Default is 0.001.

    Returns
    -------
    None
        Modifies vertices, segments, and regions lists in-place.

    Notes
    -----
    The circle is approximated by a regular polygon with n_points vertices
    distributed uniformly around the perimeter:
        theta_i = 2*pi * i / n_points,  i = 0, 1, ..., n_points-1
        x_i = center_x + radius * cos(theta_i)
        y_i = center_y + radius * sin(theta_i)

    The region marker is placed at the circle center to ensure it's inside
    the inclusion for the Triangle library.

    Typical Usage
    -------------
    Inclusions are used to model:
        - Low thermal conductivity regions (insulators)
        - High conductivity regions (heat sources)
        - Different material phases
        - Obstacles or voids in flow problems

    Examples
    --------
    Add a circular inclusion with radius 0.1m at position (0.5, 0.5):

    >>> vertices, segments, regions = [], [], []
    >>> # First create outer domain
    >>> make_box(vertices, segments, regions, 0, 0, 1, 1, attribute=1)
    >>> # Add circular inclusion
    >>> make_inclusion(vertices, segments, regions,
    ...                0.5, 0.5, 0.1, 20, attribute=100)
    >>> print(f"Total vertices: {len(vertices)}")
    Total vertices: 24

    Add multiple inclusions with different attributes:

    >>> vertices, segments, regions = [], [], []
    >>> make_box(vertices, segments, regions, 0, 0, 2, 2, attribute=1)
    >>> # High conductivity inclusion
    >>> make_inclusion(vertices, segments, regions, 0.5, 0.5, 0.2, 30,
    ...                attribute=100, max_area=0.0005)
    >>> # Low conductivity inclusion
    >>> make_inclusion(vertices, segments, regions, 1.5, 1.5, 0.15, 25,
    ...                attribute=200, max_area=0.001)

    Quality Considerations
    ----------------------
    - Use n_points >= 20 for smooth circular boundaries
    - Use smaller max_area inside inclusions for better resolution
    - Avoid overlapping inclusions (Triangle library may struggle)
    - Ensure inclusions are fully inside the outer domain

    See Also
    --------
    make_box : Add rectangular domain
    create_triangle_mesh : Complete mesh creation workflow
    """
    # Generate n_points uniformly distributed around circle
    # Note: endpoint=False ensures we don't duplicate the starting point
    theta = np.linspace(0, 2*np.pi, n_points, endpoint=False)

    # Compute x and y coordinates on circle perimeter
    xx = np.cos(theta)
    yy = np.sin(theta)

    # Get current number of vertices (for indexing new vertices)
    i = len(vertices)

    # Add circle perimeter vertices
    # Scale by radius and translate to center position
    circle_vertices = np.array([center_x + radius*xx, center_y + radius*yy]).T
    vertices.extend(circle_vertices)

    # Create segments connecting consecutive vertices around perimeter
    # Segment i connects vertex i to vertex i+1
    # Last segment wraps around: vertex (n_points-1) to vertex 0
    segment_indices = np.array([
        np.arange(i, i + n_points),         # Start vertices: [i, i+1, ..., i+n-1]
        np.arange(i + 1, i + n_points + 1)  # End vertices: [i+1, i+2, ..., i+n]
    ]).T

    # Fix last segment to close the loop (wrap around to first vertex)
    segment_indices[-1, 1] = i

    segments.extend(segment_indices)

    # Add region marker at circle center
    # Region format: [x_inside, y_inside, attribute, max_area]
    regions.append([center_x, center_y, attribute, max_area])


def create_triangle_mesh(vertices, segments, regions, mesh_options='pq33Aa'):
    """
    Generate triangle mesh using the Triangle library.

    This is a convenience wrapper around triangle.triangulate() that provides
    sensible defaults and clear documentation of the meshing process.

    Parameters
    ----------
    vertices : list
        List of [x, y] vertex coordinates defining domain boundary and features
    segments : list
        List of (i, j) tuples defining segments between vertices
    regions : list
        List of [x, y, attribute, max_area] defining regions with attributes
    mesh_options : str, optional
        Triangle library options. Default is 'pq33Aa'.
        Options:
            'p' - Triangulate a Planar Straight Line Graph
            'q' - Quality mesh (minimum angle constraint)
            '33' - Minimum angle of 33 degrees
            'A' - Assign regional attributes
            'a' - Apply maximum area constraints from regions

    Returns
    -------
    dict
        Dictionary with mesh data:
            'vertices' : ndarray, shape (nnod, 2)
                Node coordinates
            'triangles' : ndarray, shape (nel, 3)
                Element connectivity (node indices)
            'triangle_attributes' : ndarray, shape (nel,)
                Element attributes (material IDs)

    Examples
    --------
    >>> import triangle as tr
    >>> vertices, segments, regions = [], [], []
    >>> make_box(vertices, segments, regions, 0, 0, 1, 1, attribute=1)
    >>> A = dict(vertices=vertices, segments=segments, regions=regions)
    >>> B = tr.triangulate(A, 'pq33Aa')
    >>> GCOORD = B['vertices']
    >>> EL2NOD = B['triangles']
    >>> print(f"Generated {len(GCOORD)} nodes, {len(EL2NOD)} elements")

    See Also
    --------
    triangle.triangulate : Underlying triangulation function
    """
    import triangle as tr

    # Package input for Triangle library
    mesh_input = dict(
        vertices=vertices,
        segments=segments,
        regions=regions
    )

    # Generate mesh
    mesh = tr.triangulate(mesh_input, mesh_options)

    return mesh


if __name__ == "__main__":
    """
    Demonstration of mesh generation utilities.
    """
    import triangle as tr
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection

    print("=" * 70)
    print("Mesh Generation Utilities - Demonstration")
    print("=" * 70)
    print()

    # Example 1: Simple box mesh
    print("Example 1: Creating a simple box mesh")
    print("-" * 70)
    vertices, segments, regions = [], [], []
    make_box(vertices, segments, regions, 0, 0, 1, 1, attribute=1, max_area=0.01)

    A = dict(vertices=vertices, segments=segments, regions=regions)
    B = tr.triangulate(A, 'pq33Aa')

    GCOORD = B['vertices']
    EL2NOD = B['triangles']
    print(f"Generated mesh with {len(GCOORD)} nodes and {len(EL2NOD)} elements")
    print()

    # Example 2: Box with circular inclusions
    print("Example 2: Creating mesh with circular inclusions")
    print("-" * 70)
    vertices, segments, regions = [], [], []

    # Outer domain
    make_box(vertices, segments, regions, -1, -1, 2, 2, attribute=1, max_area=0.01)

    # Add 3 circular inclusions
    make_inclusion(vertices, segments, regions, -0.5, 0, 0.2, 20, attribute=100, max_area=0.001)
    make_inclusion(vertices, segments, regions, 0.5, 0, 0.2, 20, attribute=100, max_area=0.001)
    make_inclusion(vertices, segments, regions, 0, 0.5, 0.15, 20, attribute=100, max_area=0.001)

    A = dict(vertices=vertices, segments=segments, regions=regions)
    B = tr.triangulate(A, 'pq33Aa')

    GCOORD = B['vertices']
    EL2NOD = B['triangles']
    Phases = B['triangle_attributes'].flatten()

    print(f"Generated mesh with {len(GCOORD)} nodes and {len(EL2NOD)} elements")
    print(f"Region attributes: {np.unique(Phases)}")
    print()

    # Visualization (optional - requires matplotlib)
    print("Creating visualization...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Plot 1: Simple mesh
    ax1.triplot(GCOORD[:, 0], GCOORD[:, 1], EL2NOD, 'k-', linewidth=0.5)
    ax1.set_aspect('equal')
    ax1.set_title('Mesh with Inclusions')
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')
    ax1.grid(True, alpha=0.3)

    # Plot 2: Color by region attribute
    for iel, nodes in enumerate(EL2NOD):
        coords = GCOORD[nodes]
        color = 'lightblue' if Phases[iel] == 1 else 'salmon'
        triangle = Polygon(coords, facecolor=color, edgecolor='black', linewidth=0.5)
        ax2.add_patch(triangle)

    ax2.set_xlim(GCOORD[:, 0].min()-0.1, GCOORD[:, 0].max()+0.1)
    ax2.set_ylim(GCOORD[:, 1].min()-0.1, GCOORD[:, 1].max()+0.1)
    ax2.set_aspect('equal')
    ax2.set_title('Mesh colored by Material ID')
    ax2.set_xlabel('x [m]')
    ax2.set_ylabel('y [m]')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('mesh_demo.png', dpi=150, bbox_inches='tight')
    print("Saved visualization to 'mesh_demo.png'")
    print()

    print("=" * 70)
    print("Demonstration complete!")
    print("=" * 70)
