#!/bin/sh
 cd ${0%/*} || exit 1    # Run from this directory
 
 # Source tutorial run functions
 . $WM_PROJECT_DIR/bin/tools/RunFunctions
 
 application=`getApplication`
 
 ./clean.sh
runApplication blockMesh
runApplication snappyHexMesh -overwrite
runApplication checkMesh -allTopology -allGeometry

transformPoints -scale "(1e-6 1e-6 1e-6)"

runApplication $application