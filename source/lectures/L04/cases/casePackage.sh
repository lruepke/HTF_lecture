function zipCase()
{
    caseDir=$1 
    casePath=${PWD}/`dirname "$0"`
    # 1.
    cd $casePath
    rm ${caseDir}.zip
    zip -r ${caseDir}.zip $caseDir -x ${caseDir}/*00\* ${caseDir}/jupyter/data*\* 
}

zipCase Jupp_Schultz