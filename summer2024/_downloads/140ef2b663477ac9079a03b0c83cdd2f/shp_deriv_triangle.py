import numpy as np

def shp_deriv_triangle(local_coords, num_nodes):
    """
    Derivatives of triangular shape functions.

    Args:
        local_coords (numpy.ndarray): (2,nip) Local coordinates at integration points.
        num_nodes (int): Number of nodes per element.

    Returns:
        dict: A dictionary containing:
            - 'shape_functions': Shape functions at integration points (num_integration_points, num_nodes).
            - 'shape_func_deriv': Shape function derivatives with respect to xi and eta at integration points (num_integration_points, 2, num_nodes).
    """
    # Number of integration points
    num_integration_points = local_coords.shape[0]

    # Initialize arrays for shape functions and their derivatives
    shape_functions = np.zeros((num_integration_points, num_nodes))
    shape_func_deriv = np.zeros((num_integration_points, 2, num_nodes))


    for i in range(num_integration_points):
        eta2 = local_coords[i, 0]
        eta3 = local_coords[i, 1]
        eta1 = 1 - eta2 - eta3


        if num_nodes == 3:  # Linear shape functions
            shape_functions[i, 0] = eta1
            shape_functions[i, 1] = eta2
            shape_functions[i, 2] = eta3

            shape_func_deriv[i, 0, 0] = -1
            shape_func_deriv[i, 0, 1] = 1
            shape_func_deriv[i, 0, 2] = 0

            shape_func_deriv[i, 1, 0] = -1
            shape_func_deriv[i, 1, 1] = 0
            shape_func_deriv[i, 1, 2] = 1

        elif num_nodes == 6:  # Quadratic shape functions
            shape_functions[i, 0] = eta1*(2*eta1-1)
            shape_functions[i, 1] = eta2*(2*eta2-1)
            shape_functions[i, 2] = eta3*(2*eta3-1)
            shape_functions[i, 3] = 4*eta2*eta3
            shape_functions[i, 4] = 4*eta1*eta3
            shape_functions[i, 5] = 4*eta1*eta2

            shape_func_deriv[i, 0, 0] = 1-4*eta1 -1+4*eta2
            shape_func_deriv[i, 0, 1] = 0
            shape_func_deriv[i, 0, 2] = 0
            shape_func_deriv[i, 0, 3] = 4*eta3
            shape_func_deriv[i, 0, 4] = -4*eta3
            shape_func_deriv[i, 0, 5] = 4*eta1-4*eta2

            shape_func_deriv[i, 1, 0] = 1-4*eta1
            shape_func_deriv[i, 1, 1] = 0
            shape_func_deriv[i, 1, 2] = -1+4*eta3
            shape_func_deriv[i, 1, 3] = 4*eta2
            shape_func_deriv[i, 1, 4] =4*eta1-4*eta
            shape_func_deriv[i, 1, 5] = -4*eta2

        elif num_nodes == 7:  # Quadratic shape functions with internal node
            shape_functions[i, 0] = eta1*(2*eta1-1)+ 3*eta1*eta2*eta3
            shape_functions[i, 1] = eta2*(2*eta2-1)+ 3*eta1*eta2*eta3
            shape_functions[i, 2] = eta3*(2*eta3-1)+ 3*eta1*eta2*eta3
            shape_functions[i, 3] = 4*eta2*eta3 - 12*eta1*eta2*eta3
            shape_functions[i, 4] = 4*eta1*eta3 - 12*eta1*eta2*eta3
            shape_functions[i, 5] = 4*eta1*eta2 - 12*eta1*eta2*eta3
            shape_functions[i, 6] = 27*eta1*eta2*eta3

            shape_func_deriv[i, 0, 0] = 1-4*eta1+3*eta1*eta3-3*eta2*eta3
            shape_func_deriv[i, 0, 1] = -1+4*eta2+3*eta1*eta3-3*eta2*eta3
            shape_func_deriv[i, 0, 2] = 3*eta1*eta3-3*eta2*eta3
            shape_func_deriv[i, 0, 3] = 4*eta3+12*eta2*eta3-12*eta1*eta3
            shape_func_deriv[i, 0, 4] = -4*eta3+12*eta2*eta3-12*eta1*eta3
            shape_func_deriv[i, 0, 5] = 4*eta1-4*eta2+12*eta2*eta3-12*eta1*eta3
            shape_func_deriv[i, 0, 6] = -27*eta2*eta3+27*eta1*eta3

            shape_func_deriv[i, 1, 0] = 1-4*eta1+3*eta1*eta2-3*eta2*eta3
            shape_func_deriv[i, 1, 1] = +3*eta1*eta2-3*eta2*eta3
            shape_func_deriv[i, 1, 2] = -1+4*eta3+3*eta1*eta2-3*eta2*eta3
            shape_func_deriv[i, 1, 3] = 4*eta2-12*eta1*eta2+12*eta2*eta3
            shape_func_deriv[i, 1, 4] = 4*eta1-4*eta3-12*eta1*eta2+12*eta2*eta3
            shape_func_deriv[i, 1, 5] = -4*eta2-12*eta1*eta2+12*eta2*eta3
            shape_func_deriv[i, 1, 6] = 27*eta1*eta2-27*eta2*eta3

    # Package the results in a dictionary
    results = {
        'shape_functions': shape_functions,
        'shape_func_deriv': shape_func_deriv
    }

    return results

# Example usage
#local_coords = np.array([[0.333, 0.2, 0.6], [0.333, 0.6, 0.2]])
#num_nodes = 3  # for a linear triangular element
#results = shp_deriv_triangle(local_coords, num_nodes)

# Get the shape functions and derivatives for the first integration point
#integration_point_index = 0
#shape_funcs_at_point = results['shape_functions'][integration_point_index]
#shape_func_derivs_at_point = results['shape_func_deriv'][integration_point_index]

#print("Shape functions at integration point:")
#print(shape_funcs_at_point)
#print("\nShape function derivatives at integration point:")
#print(shape_func_derivs_at_point)
