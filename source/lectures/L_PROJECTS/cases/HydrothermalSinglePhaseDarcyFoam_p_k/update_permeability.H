// 1. update internal value
// Iyer, K., Schmid, D. W., Planke, S., & Millett, J. (2017). Modelling hydrothermal venting in volcanic sedimentary basins: Impact on hydrocarbon maturation and paleoclimate. Earth and Planetary Science Letters, 467, 30–42. https://doi.org/10.1016/j.epsl.2017.03.023
{
    forAll(p,i)
    {
        if(p[i]>pl[i])
        {
            permeability[i]=permeability[i]*pow(p[i]/pl[i], 2); //Eq.(7) of Iyer et al.(2017)
            if(permeability[i]>1E-13)
            {
                permeability[i] = 1E-13;
            }
        }
        else
        {
            permeability[i]=permeability0[i];
        }
        
    }
    forAll(mesh.cellZones(), zoneI)
        {
            const cellZone& cz = mesh.cellZones()[zoneI];
            const labelList& cells = cz;

            if (cz.name() == "intrusion")
            {
                forAll(cells, iCells)
                {   
                    const label cellId = cells[iCells];
                    permeability[cellId]=permeability0[cellId];
                }
            }
        }




}