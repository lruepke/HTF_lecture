// 1. update internal value
// Aarnes, I., Svensen, H., Connolly, J. A. D., & Podladchikov, Y. Y. (2010). How contact metamorphism can trigger global climate changes: Modeling gas generation around igneous sills in sedimentary basins. Geochimica et Cosmochimica Acta, 74(24), 7179–7195. https://doi.org/10.1016/j.gca.2010.09.011
{
    double Lc = 320; //Latent heat of crystallization, [kJ/kg]
    double Tl = 1423; //Liquidus temperature, [K]
    double Ts = 1223; //Solidus temperature, [K]
    double deltaT_LS = Tl-Ts; 
    double cp0 = cp0_rock.value();
    forAll(T,i)
    {
        if(T[i]>Ts && T[i]<Tl)
        {
            cp_rock[i] = cp0*(1 + Lc/(deltaT_LS*cp0)); //eq(2) of Arnes et al.(2010)
        }else
        {
            cp_rock[i] = cp0;
        }
    }
}