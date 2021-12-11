Cl_print=zeros(nl_soil,365*how_many_year);
Cb_print=zeros(nl_soil,365*how_many_year);
Ch_print=zeros(nl_soil,365*how_many_year);
Nl_print=zeros(nl_soil,365*how_many_year);
CNl_print=zeros(nl_soil,365*how_many_year);

Amm_print=zeros(nl_soil,365*how_many_year);
Nit_print=zeros(nl_soil,365*how_many_year);

Drain_Nit_print=zeros(1,365*how_many_year);
Drain_Amm_print=zeros(1,365*how_many_year);
    
if from_previous_data ~= 0
    UP_amm_print=zeros(nl_soil,365*how_many_year);
    UP_amm_all_print=zeros(nl_soil,365*how_many_year);
    UP_amm_all_m2_print=zeros(nl_soil,365*how_many_year);    
    UP_nit_print=zeros(nl_soil,365*how_many_year);
    UP_nit_all_print=zeros(nl_soil,365*how_many_year);
    UP_nit_all_m2_print=zeros(nl_soil,365*how_many_year);    
    UP_N_m2_print=zeros(1,365*how_many_year);
    UP_nit_m2_print=zeros(nl_soil,365*how_many_year);
    UP_amm_m2_print=zeros(nl_soil,365*how_many_year);    
    UP_amm_active_print=zeros(nl_soil,365*how_many_year);
    UP_nit_active_print=zeros(nl_soil,365*how_many_year);

    LCH_amm_print=zeros(nl_soil,365*how_many_year);
    LCH_nit_print=zeros(nl_soil,365*how_many_year);
    LCH_N_m2_print=zeros(1,365*how_many_year);    
    TLCH_amm_m2_print=zeros(1,365*how_many_year);
    TLCH_nit_m2_print=zeros(1,365*how_many_year);
    LCH_nit_m2_print=zeros(nl_soil,365*how_many_year);
    LCH_amm_m2_print=zeros(nl_soil,365*how_many_year);

    MIN_net_print=zeros(nl_soil,365*how_many_year);    
    MIN_gross_print=zeros(nl_soil,365*how_many_year);
    IMM_net_print=zeros(nl_soil,365*how_many_year); 
    IMM_gross_print=zeros(nl_soil,365*how_many_year);
    Nreg_print=zeros(nl_soil,365*how_many_year);
    DECl_print=zeros(nl_soil,365*how_many_year);
    PHI_print=zeros(nl_soil,365*how_many_year);
    phi_print=zeros(nl_soil,365*how_many_year); 
    fSd_print=zeros(nl_soil,365*how_many_year);
    fTd_print=zeros(nl_soil,365*how_many_year);
    mberrorN_print=zeros(1,365*how_many_year);
    mberrorC_print=zeros(1,365*how_many_year);

    flow_Nit_m2_print=zeros(nl_soil-1,365*how_many_year);
    flow_Amm_m2_print=zeros(nl_soil-1,365*how_many_year);

    if SWITCHES.CN.Bioturbation
        Cbiorate_print=zeros(1,365*how_many_year);    
        Nbiorate_print=zeros(1,365*how_many_year);   
        bioNerror_print=zeros(1,365*how_many_year);
        bioCerror_print=zeros(1,365*how_many_year);

        Cl_bio_change_print=zeros(nl_soil,365*how_many_year);                
        Cin_m3_print=zeros(nl_soil,365*how_many_year);                                            
        Cout_m3_print=zeros(nl_soil,365*how_many_year); 
    end

    if SWITCHES.Active_Nuptake == 1;
        Amm_DEM_print=zeros(nl_soil,365*how_many_year);
        Nit_DEM_print=zeros(nl_soil,365*how_many_year);
    end
end