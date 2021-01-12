#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

application=`getApplication`

./clean.sh
# # option 1 regular box
# runApplication blockMesh
# option 2 gmsh unstructured mesh
gmsh gmsh/mesh.geo -3 -o gmsh/mesh.msh -format msh22
gmshToFoam gmsh/mesh.msh
changeDictionary
# try to comment renumberMesh and see the matrix 
renumberMesh -overwrite

# runApplication $application
# foamToVTK -useTimeName