function zipCase()
{
    caseDir=$1 
    casePath=${PWD}/`dirname "$0"`
    # 1.
    cd $casePath
    rm ${caseDir}.zip
    zip -r ${caseDir}.zip $caseDir -x ${caseDir}/*00\* ${caseDir}/jupyter/data*\* ${caseDir}/ventT.* ${caseDir}/VTK\*
}
zipCase fault_flow
zipCase cooling_intrusion
zipCase HydrothermalSinglePhaseDarcyFoam_Cpr
zipCase HydrothermalSinglePhaseDarcyFoam_p_k
zipCase sill_instrusion