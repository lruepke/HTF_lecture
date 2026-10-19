
function writelines
{
    fname_data=$1
    echo "pts_${fname_data}={}; lines_${fname_data}={};" >${fname_data}.geo
    awk -v name=$fname_data 'NR<2{print "pts_"name"[#pts_"name"[]]=newp;Point(pts_"name"[#pts_"name"[]-1])={"$1*1000","$2*1000",zmin,lc};"}' ${fname_data}.dat >>${fname_data}.geo
    awk -v name=$fname_data 'NR>1{print "pts_"name"[#pts_"name"[]]=newp;Point(pts_"name"[#pts_"name"[]-1])={"$1*1000","$2*1000",zmin,lc}; lines_"name"[#lines_"name"[]]=newl; Line(lines_"name"[#lines_"name"[]-1])={pts_"name"[#pts_"name"[]-2], pts_"name"[#pts_"name"[]-1]};"}' ${fname_data}.dat >>${fname_data}.geo
}

writelines Tertiary
writelines UpperCretaceous
writelines sill1
writelines sill2
writelines sill3
writelines sill4
writelines sill5