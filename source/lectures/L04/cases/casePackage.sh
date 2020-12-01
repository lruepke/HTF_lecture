caseDir=Jupp_Schultz 
casePath=${PWD}/`dirname "$0"`
rm ${casePath}/${caseDir}.zip
zip -r ${casePath}/${caseDir}.zip ${casePath}/$caseDir -x ${casePath}/${caseDir}/*00\* ${casePath}/${caseDir}/jupyter/data*\* 