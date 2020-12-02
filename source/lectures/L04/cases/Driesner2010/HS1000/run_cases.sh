function run_group()
{
    group=$1
    for k_df in 1E-16 5E-16 1E-15 2E-15 4E-15 6E-15 8E-15 1E-14 2E-14 4E-14 6E-14 8E-14 1E-13
    do 
        caseDir=DF_${k_df}
        python run_case.py $group/$caseDir
    done
}

run_group Pipe_1E-14
run_group Pipe_2E-14
run_group Pipe_1E-16
run_group Pipe_5E-15