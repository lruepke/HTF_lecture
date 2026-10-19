#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

application=`getApplication`

./clean.sh
# 1. meshing using gmsh
gmsh gmsh/make_mesh.geo -3 -o gmsh/intrusion_mesh.msh -format msh22

# 2.convert gmsh (.msh file) to openfoam format (polyMesh)
gmshToFoam gmsh/intrusion_mesh.msh

# 3..change front and back patches as empty to form a 2D mesh
changeDictionary

# 4. set permeability
runApplication setFields
# 5. run
runApplication $application
# 6. postprocessing
# ./postProcess.sh