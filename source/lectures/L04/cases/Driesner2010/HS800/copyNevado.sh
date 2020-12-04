baseDir_Nevado=/Data/zguo/HTF_lecture/source/lectures/L04/cases/Driesner2010/HS1000

for k_df in 1E-16 3E-16 5E-16 8E-16 1E-15 2E-15 3E-15 4E-15 6E-15 8E-15 1E-14 3E-14 6E-14 8E-14 1E-13 2E-13 3E-13 4E-13 1E-12
do 
    scp nevado:${baseDir_Nevado}/${k_df}/ventT.txt results/ventT_${k_df}.txt
done
