#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

application=`getApplication`

# 1. clean old results
# ./clean.sh
# 2. mesh
blockMesh
# 5. run 
runApplication $application