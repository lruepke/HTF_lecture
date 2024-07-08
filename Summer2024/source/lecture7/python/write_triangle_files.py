def write_mesh_files(mesh_data, model_name):
    
    node_filename = f"{model_name}.node"
    ele_filename = f"{model_name}.ele"
    
    # Extract GCOORD, EL2NOD, Phases, and Node_ids data
    gcoord = mesh_data.get('GCOORD')
    el2nod = mesh_data.get('EL2NOD') + 1  # Convert to 1-based indexing
    phases = mesh_data.get('Phases') + 1  # Convert to 1-based indexing
    node_ids = mesh_data.get('Node_ids')

    if gcoord is None:
        raise ValueError("The mesh data does not contain 'GCOORD' key.")
    if el2nod is None:
        raise ValueError("The mesh data does not contain 'EL2NOD' key.")
    if phases is None:
        raise ValueError("The mesh data does not contain 'Phases' key.")
    if node_ids is None:
        raise ValueError("The mesh data does not contain 'Node_ids' key.")

    # Calculate the number of nodes to omit from the end of GCOORD
    num_nodes_to_omit = el2nod.shape[0]

    # Adjust the GCOORD and Node_ids to omit the last num_nodes_to_omit entries
    gcoord = gcoord[:-num_nodes_to_omit]
    node_ids = node_ids[:-num_nodes_to_omit]

    # Write the .node file
    with open(node_filename, 'w') as node_file:
        # Write the header line
        node_file.write(f"{len(gcoord)} 2 0 1\n")

        # Write each point with the corresponding id from Node_ids
        for index, (coords, node_id) in enumerate(zip(gcoord, node_ids), start=1):
            x = coords[0]
            y = coords[1]
            marker = node_id
            node_file.write(f"{index} {x} {y} {marker}\n")

    # Write the .ele file
    with open(ele_filename, 'w') as ele_file:
        num_elements = len(el2nod)
        num_nodes_per_element = min(len(el2nod[0]), 6)  # Use up to 6 nodes per element

        # Write the header line
        ele_file.write(f"{num_elements} {num_nodes_per_element} 1\n")

        # Write each element
        for i, element in enumerate(el2nod, start=1):
            nodes = " ".join(map(str, element[:6]))  # Use only the first 6 nodes
            phase_id = phases[i - 1]
            ele_file.write(f"{i} {nodes} {phase_id}\n")