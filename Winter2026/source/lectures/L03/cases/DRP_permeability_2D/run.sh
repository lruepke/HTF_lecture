#!/bin/sh
 cd ${0%/*} || exit 1    # Run from this directory
 
 # Source tutorial run functions
 . $WM_PROJECT_DIR/bin/tools/RunFunctions
 
 application=`getApplication`
 
 ./clean.sh
runApplication blockMesh
runApplication snappyHexMesh -overwrite
# 2-D case: snappyHexMesh needs a 3-D mesh, so frontAndBack is a patch until here
foamDictionary constant/polyMesh/boundary -entry entry0/frontAndBack/type -set empty > /dev/null
runApplication checkMesh -allTopology -allGeometry

transformPoints "scale=(1e-6 1e-6 1e-6)"

runApplication $application