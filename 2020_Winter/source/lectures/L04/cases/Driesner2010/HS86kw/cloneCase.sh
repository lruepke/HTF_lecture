
baseGroup=1E-16

for k_df in 3E-16 5E-16 
do 
    newCaseDir=${k_df}
    rm -rf ${newCaseDir}
    cp -rf ${baseGroup} ${newCaseDir}
    awk -v perm=${k_df}  '{gsub(/internalField   uniform 1e-16;/, "internalField   uniform "perm";"); print }' ${newCaseDir}/0/permeability > permeability.tmp
    mv permeability.tmp ${newCaseDir}/0/permeability
done

for k_df in 8E-16 1E-15 2E-15 3E-15 4E-15 6E-15 8E-15 1E-14 3E-14 6E-14 8E-14 1E-13 2E-13 3E-13 4E-13 1E-12
do 
    newCaseDir=${k_df}
    rm -rf ${newCaseDir}
    cp -rf ${baseGroup} ${newCaseDir}
    awk -v perm=${k_df}  '{gsub(/internalField   uniform 1e-16;/, "internalField   uniform "perm";"); print }' ${newCaseDir}/0/permeability > permeability.tmp
    mv permeability.tmp ${newCaseDir}/0/permeability
    # change endTime
    awk '{gsub(/endTime 757680000000;/, "endTime 257680000000;"); print }' ${newCaseDir}/system/controlDict > controlDict.tmp
    mv controlDict.tmp ${newCaseDir}/system/controlDict
done
