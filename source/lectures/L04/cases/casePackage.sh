function zipCase()
{
    caseDir=$1 
    casePath=${PWD}/`dirname "$0"`
    # 1.
    cd $casePath
    rm ${caseDir}.zip
    zip -r ${caseDir}.zip $caseDir -x ${caseDir}/*00\* ${caseDir}/jupyter/data*\* ${caseDir}/ventT.* ${caseDir}/VTK\*
}
function zipCaseDriesner()
{
    caseDir=$1 
    casePath=${PWD}/`dirname "$0"`
    # 1.
    cd $casePath
    rm ${caseDir}.zip
    zip -r ${caseDir}.zip $caseDir 
}
zipCase Jupp_Schultz
zipCaseDriesner Driesner2010