#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

application=`getApplication`

# 1. clean old results
./clean.sh
# 2. generate mesh using gmsh
# gmsh gmsh/mesh.geo -3 -o gmsh/mesh.msh -format msh22
# # 3. convert gmsh format to OpenFOAM format
# gmshToFoam gmsh/mesh.msh 
# # 4. set empty patch for 2D case
# changeDictionary
blockMesh
# 5. run 
runApplication $application