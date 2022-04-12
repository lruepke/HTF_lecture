#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

application=`getApplication`

# 1. slice 
postProcess -func surfaces

# # 2. profile sampling
postProcess -func sampleDict -latestTime

# 3. mass flux profile on patch
# $application -postProcess -dict system/massFluxDict -latestTime

# 4. last time
# foamToVTK -noInternal -useTimeName -latestTime

# 5. temperature on seafloor patch
# foamToVTK -noInternal -useTimeName  -fields '(T)' -excludePatches '(frontAndBack heatsource bottom sidewalls)'