%%-----------------------------------------------------------------------%%
%%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>%%
%+             Calculate Dynamics of Soil Carbon and Nitrogen            +%
%%-----------------------------------------------------------------------%%
%%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>%%
%
%   This model is used to run Soil Carbon and Nitrogen cycle
%   presented by Porpoarato et al (2003).
%
%%-----------------------------------------------------------------------%%
%
%   Created by      :   Juan Quijano (mainframe) & Dongkook Woo (nitrogen
%   age & agricultural part) & Susana Roque-Malo (exudation, rhizosphere,
%   crunch coupling part)
%   Written in      :   Matlab
%   Date            :   10-02-14
%
%%-----------------------------------------------------------------------%%
% This code is subjected to Corn-corn-soybean rotation
%%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>%%

clc
clear all
close all
tic

CORN = 1; % 0 if soybean
% load basic parameters
if CORN == 1
    load cndatacorn.mat;
    start_year=2008;%1908;%%2006;%2008;%
    end_year=2008;%2107;%2007;%2012;%2007;%2011;%%2050;%1909;%2047%2011%
else
    load cndatasoy.mat;
    start_year=2010;%1908;%%2006;%2008;%
    end_year=2010;%2107;%2007;%2012;%2007;%2011;%%2050;%1909;%2047%2011%

end

% %-----------------------------------------------------------------------%%
% Setting before running the model                                        %
% %-----------------------------------------------------------------------%%
% 1. Initialization with corn/soy = 1
initialization=0;
% 
% 2. From the VERY 1st run, choose 0, otherwise, choose 1
from_previous_data =0;

    if from_previous_data==1
        load('COSOY_initial_run_new.mat','VARIABLES','PARAMS') % Min & Imm are coupled
        
        Previous_Cl=VARIABLES.Cl;
        Previous_Cb=VARIABLES.Cb;
        Previous_Ch=VARIABLES.Ch;
        Previous_Amm=VARIABLES.Amm;
            Amm=VARIABLES.Amm;

        Previous_Amm_mo=VARIABLES.Amm_mo;
            Amm_mo=VARIABLES.Amm_mo;
        Previous_Amm_im=VARIABLES.Amm_im;
            Amm_im=VARIABLES.Amm_im;
            
            
        Previous_NH3=VARIABLES.NH3;
            NH3=VARIABLES.NH3;
        Previous_Nit=VARIABLES.Nit;
            Nit=VARIABLES.Nit;
        Previous_CNl=VARIABLES.CNl;
        Previous_Nl=VARIABLES.Nl;

        Previous_kl=PARAMS.CN.kl;
        Previous_kd=PARAMS.CN.kd;
        Previous_kh=PARAMS.CN.kh;
    end
% 
% 3. want to save "save_whole_cycle" = 1 otherwise =0
save_whole_cycle=1;
% 
% 4. fixed kl, kd, kh
% For fixed kl, kd, and kh, choose 1, otherwise, choose 0
% Caution1: Before using 'fixed_kl_kd_kh', first you need to run the model with
% the given years
% Caution2: from_previous_data = 1, it is not necessary to put number here
fixed_kl_kd_kh = 0;
% 
    if fixed_kl_kd_kh==1
        load('kl_kd_kh.mat','kl','kd','kh'), % save('kl_kd_kh','kl','kd','kh')
        Fixed_kl=kl;
        Fixed_kd=kd;
        Fixed_kh=kh;
    end
%         
% 6. Fix biomass drops & crop carbon & crop water uptake = 1; otherwise =0
Fix_Crop_Vaues=0;     

    if Fix_Crop_Vaues==1
        load('Corn_S1.mat','TBMLa','TBMLaNoLitterBack','RootDeathWhenHarvest','TBMLaOnlyLitterBack','TBMLaNoLitterBack',...
            'CARBONS',...
            'adupsoil_layer_day_point')
        
        pre_TBMLa=TBMLa;
        pre_TBMLaNoLitterBack=TBMLaNoLitterBack;
        pre_RootDeathWhenHarvest=RootDeathWhenHarvest;
        pre_TBMLaOnlyLitterBack=TBMLaOnlyLitterBack; 
        pre_TBMLaNoLitterBack=TBMLaNoLitterBack;
        
        pre_grain_carbon_day_point=CARBONS.grain_carbon_day_point;
        pre_leaf_carbon_day_point=CARBONS.leaf_carbon_day_point;
        pre_rhizome_carbon_day_point=CARBONS.rhizome_carbon_day_point;
        pre_root_carbon_day_point=CARBONS.root_carbon_day_point;
        pre_stem_carbon_day_point=CARBONS.stem_carbon_day_point;
        
        pre_adupsoil_layer_day_point=adupsoil_layer_day_point;
    end
    
% 5. Choose scenarios 1 to 4 for ppt change, 5 for 100 year simulation
scenarios=1;
% 
% 6. Added Carbon change factor [g/m2] increase / FOR calibration!! /
% NP:1, Corn:1, SG:1, MG:1
% This is always "1" if you do not want to modify the carbon change!
% This value will be multiplied to 'plant total carbon change'
Add_Carbon=1;
% 
% 7. Choose the start, and end years, and soybean years!
% and first day of corn and soy day.
% start_year=2010;%1908;%%2006;%2008;%
% end_year=2010;%2107;%2007;%2012;%2007;%2011;%%2050;%1909;%2047%2011%
    PreSoy=1809:2:2007;
    PostSoy=2010:3:2107;
    PreCorn=1808:2:2007;
    PostCorn1=2008:3:2107;
    PostCorn2=2009:3:2107;
soy_year=[PreSoy PostSoy];
corn_year=[PreCorn PostCorn1 PostCorn2];
% 
% 8. How many years do you want to run this model? / total 43 years *24->maximum in a run. 
how_many_year=end_year-start_year+1;
%       
% 9. Choose which crops? 1 for NP, 2 for Corn/Soy, 3 for SG, and 4 for MG
choose_crop=2;
% 
% %-----------------------------------------------------------------------%%
% Tuning parameter                                                        %
% %-----------------------------------------------------------------------%%
% 1. Only kn: rate of nitrification [m3/d/gN] / 0.003;
% Corn:0.00004, SG:0.0002, MG: 0.00021,
PARAMS.CN.kn=0.00004;%
kn=PARAMS.CN.kn;

% 2. N uptake parameters. F: rescaled diffusion coefficient
% , 0.1 [m/m3] (D'Ddoricoe et al., 2003). d: the nonlinear dependence of
% the diffusion coefficient, 3 (Porporato et al., 2003).
% Corn:F=0.06, SG:0.015, MG: 0.035
PARAMS.F=0.06;% merely chosen
PARAMS.d=3;% Typical Value
% 
% 3. rh: the decomposed litter that undergoes humification,
% which is in the range 0.15-0.35 depending on the composition of plant
% residues.
% Corn:0.15, SG: 0.15, MG:0.15
PARAMS.CN.rhmin=0.15;%0.15
% 
% 4. a site specific parameter that needs to be estimated with observed
% N2O data. from decomposition*K2. mean values is 0.02
% Dnitrification ratio is based on average several average data, thus,
% reduce of increase based on a fraction. and the paper said that these
% values should be well matched up with observed values.
% Corn:0.02, SG: 0.02, MG:0.02
PARAMS.K2=0.02;
% 
% 5. N remobilization: should be always 1. You can control this with
% 'Nremobilization' parameters
SWITCHES.Nremobil=1; %SUSANA EDIT - MICROBIAL ACTIVITY HANDLED BY CRUNCH
% 
% %-----------------------------------------------------------------------%%
% From published papers 
% %-----------------------------------------------------------------------%%
% 1. Different parameters and varibles for each crop
% Caution: Choose the plant functional type. 13 for Soy, 14 for Corn, 
%      15 for NP, 16 for MG, and 17 for SG     
% Caution: Choose the points? 
% 16 by 16 > 68 for NP, 204 for Corn/Soy, 76 for SG, and 196 for MG
%  8 by  8 > 19 for NP,  46 for Corn/Soy, 22 for SG, and  43 for MG
% 
% 2. Corn-Corn-Soy rotation
if choose_crop==2
    choose_point=46;
    choose_pft_Corn=14;
    choose_pft_Soy=13;
%     1. root fraction [%] 
%     see below

%     2. root death/litterfall [-]
    CR_ratio_Corn=0.21;
    CR_ratio_Soy=0.35;
    PARAMS.CR_ratio_Corn=CR_ratio_Corn;
    PARAMS.CR_ratio_Soy=CR_ratio_Soy;
%     
%     3. leaf life span [day]
    leafspan_Corn=115+22; 
    leafspan_Soy=65; 
%     
%     4. Abovegound biomass back to soil when harvests [%]
%     Percent of aboveground C (stem + leaf! no grain)
    HavestResidueRate_Corn=0.5;
    HavestResidueRate_Soy=1;
%     
%     5. C/N ratio for aboveground and belowground 
%     Caution: residue-> corn: 60-100, and soybean residue : 20 - 40
%     But from EBI site, above corn 30.7-32.2, soy, 9
%     Root CN
    CN_leaf_Corn=22;
    CN_stem_Corn=56;
    CN_grain_Corn=CN_stem_Corn;
    CN_root_Corn=27;
%     
    CN_leaf_Soy=6; 
    CN_stem_Soy=22;
    CN_grain_Soy=CN_stem_Soy;
    CN_root_Soy=13;
%     
%     6. N fixation rate = N uptake from N fixation / N uptake from soil
%     corn-corn-soybean
    nfix=[0.00; 0.00; 0.3];
%     
%     8. N remobilization
%     08 09 10 11 lest
    Nremobilization=[0;0;0;0;0];
end
% 


% %-----------------------------------------------------------------------%%
% From study site                                                         %
% %-----------------------------------------------------------------------%%
% 1. defalt values
% 1.1. Fertilizer [lb/acre]. Caution: check 3 years from the start year
% Only for Corn and Switchgrass
% Caution: below values will be changed due to the start day. if Corn Soy Corn.
% Caution: Do not put "0" for "Day_1,2, and 3rd", instead just put "1"
% 6 years, 3 year of observation, 3 years of future senarios
Nitrate_fertilizer_amount_1st_year=0;
Nitrate_fertilizer_amount_2nd_year=0;
Nitrate_fertilizer_amount_3rd_year=0;
Nitrate_fertilizer_amount_4th_year=0;
Nitrate_fertilizer_amount_5th_year=0;
Nitrate_fertilizer_amount_6th_year=0;

Ammonium_fertilizer_amount_1st_year=0;
Ammonium_fertilizer_amount_2nd_year=0;
Ammonium_fertilizer_amount_3rd_year=0;
Ammonium_fertilizer_amount_4th_year=0;
Ammonium_fertilizer_amount_5th_year=0;
Ammonium_fertilizer_amount_6th_year=0;

Organic_N_fertilizer_1st_year=0;
Organic_N_fertilizer_2nd_year=0;
Organic_N_fertilizer_3rd_year=0;
Organic_N_fertilizer_4th_year=0;
Organic_N_fertilizer_5th_year=0;
Organic_N_fertilizer_6th_year=0;

Day_of_fertilzer_1st_year=1;
Day_of_fertilzer_2nd_year=1;
Day_of_fertilzer_3rd_year=1;%
Day_of_fertilzer_4th_year=1;
Day_of_fertilzer_5th_year=1;
Day_of_fertilzer_6th_year=1;
% 
% 2. Atmospheric N deposition (mean value over 2002 - 2011)
% 2008 2009 2010 2011 and then average
Deposition_NH4_N=[0.31/365; 0.36/365; 0.22/365; 0.28/365; 0.29/365]; % g/m2/day
Deposition_NO3_N=[0.25/365; 0.25/365; 0.19/365; 0.19/365; 0.24/365]; % g/m2/day

% 
% 3. Set number of soil layers (from MLCan simulation)
PARAMS.nl_soil=12;
nl_soil=PARAMS.nl_soil;
N = nl_soil;
% 
% 4. Soil pH
VARIABLES.pH=6.*ones(nl_soil,1); % Smith et al., 2012
% 
% 4. Different parameters and varibles for each crop
% 4.1 Corn-Corn-Soy rotation
if choose_crop==2
%     SUSANA EDIT -REMOVE PALMS
% soil type, classifications in Woo 2017
    soiltype=[2015 2015 2015 2015 933,...
        933 933 3234 3234 4218,...
        4218 4218];
% 
%     Nitrogen fertilizer amount [g/m2], UAN, 25% NH3-, 25% NH4+, 50% Urea
%     Volatilization of UAN, % OF N applied 14-37 => average 25.5%
%     Volatilization of Urea, % OF N applied 10-31 => average 20.5%
%     Thus NH4+'s vlolatilization is 0.61%
%     However, from now on, Volatilization model is developed as well as urea model.
%     Therefore, we do not need to assume the volatilization from
%     fertilizer applicaton.
    Nitrate_fertilizer_amount_1st_year=(16.8*0.25);
    Nitrate_fertilizer_amount_2nd_year=(20.2*0.25);
    Nitrate_fertilizer_amount_3rd_year=0;%
    Nitrate_fertilizer_amount_4th_year=(18.0*0.25);
    Nitrate_fertilizer_amount_5th_year=(18.0*0.25);
    Nitrate_fertilizer_amount_6th_year=0;
    
    Ammonium_fertilizer_amount_1st_year=(16.8*0.25);
    Ammonium_fertilizer_amount_2nd_year=(20.2*0.25);
    Ammonium_fertilizer_amount_3rd_year=0;%
    Ammonium_fertilizer_amount_4th_year=(18.0*0.25);
    Ammonium_fertilizer_amount_5th_year=(18.0*0.25);
    Ammonium_fertilizer_amount_6th_year=0;
                
    Organic_N_fertilizer_1st_year=(16.8*0.5);
    Organic_N_fertilizer_2nd_year=(20.2*0.5);
    Organic_N_fertilizer_3rd_year=0;
    Organic_N_fertilizer_4th_year=(18.0*0.5);
    Organic_N_fertilizer_5th_year=(18.0*0.5);
    Organic_N_fertilizer_6th_year=0;
    
    Day_of_fertilzer_1st_year=126; % May 6
    Day_of_fertilzer_2nd_year=132; % May 12
    Day_of_fertilzer_3rd_year=1;
    Day_of_fertilzer_4th_year=131; % May 11
    Day_of_fertilzer_5th_year=131; % May 11
    Day_of_fertilzer_6th_year=1;
end
% 
% %
% %-----------------------------------------------------------------------%%
% Shared condition and parameters                                         %
% %-----------------------------------------------------------------------%%
% 1. If you want to use humus pool, choose '1'. Otherwise, choose '0'.
% Caution: I am not sure, whether it works without 'humus pool'
SWITCHES.CN_type = 1;
% 
% 2. If you want to use active N uptake, choose '1', 
% or to only consider pasive N uptake '0'.
% Caution: must use 'active N uptake' for sure
SWITCHES.Active_Nuptake = 1;
% 
% 3. Turn on Bioturbation=1, Turn off =0
% Caution: I am not sure, whether it works without 'bioturbation'
SWITCHES.CN.Bioturbation=1;
% 
% 4. Turn on New Bioturbation=1, Turn off =0
% Caution: If you do not want to use bioturbation at all, you need to turn
% 'SWITCHES.CN.Bioturation' off. This is the subversion of Juan's version.
% Maily because bioturbation should happen all of the soil carbon. Not only
% for the carbon in litter pool
% This only works for 3 pools!!!!
SWITCHES.CN.Bioturbation_new=1;
% 
% 4. Turn on Lignin effect on liter decomposition=1
% or, Turn off =0
% Caution: This lignin consideration works only under bioturbation. if you
%   want to improve this matter. Feel free!
SWITCHES.CN.Lignin=0;
%     
% 5. C/N ratio in biomass pool
CN_biomass=11.5;
% 
% #6 -removed
% 
% 7. Accodring to nitrate, and ammonium soluability, and charge(+ or -)
% and added ammonia, which has high affinity of water.
PARAMS.CN.a_Nit=1.00; 
PARAMS.CN.a_Amm= 0.15;% 0.05; %susana edit 10/18/2020
PARAMS.CN.a_NH3=1.00; 
% 
% 8. Because of adapting this model from Quijano.
CONSTANTS.dtime=86400;
PARAMS.CanStruc.nspecies=1;
nspecies = 1;
% 
% 10. Soil depth of each layer [mm] (from MLCan simulation)
VERTSTRUC.dzsmm=[0.06;0.08;0.09;0.1;0.1;0.17;0.2;0.2;0.4;0.4;0.45;1.25]*1000;
dz = VERTSTRUC.dzsmm./1000;%thickness of each layer in m
% 
%  
% 
% 
% 11. rr: the constnat rr defines the fraction of decomposed
% organic carbon that goes into respiration (CO2 production), which is
% usually estimated in the 0.6-0.8 interval.
PARAMS.CN.rr=0.6;
% 
% 12. Turn on Nitrogen fixation, Turn on == 1 or, Turn off =0
% Caution: This only works for soybean and miscanthus
SWITCHES.CN.Nfix=1;
% 
% 13. Turn on linear N deposition = 1, Turn off = 0. 
% Turn on == 1 or, Turn off =0
SWITCHES.CN.Ndeposition=1;
% 
% 14. Turn on Modified eqn of immobilization = 1, Turn off = 0. 
% Turn on == 1 or, Turn off =0
% From now on, I am not sure whether it works if you choose 1
SWITCHES.ModifiedEq=0;
% 
% 15. Turn on linear CN ratio change for SG, and MG. 
% Turn on == 1 or, Turn off =0
% From now on, I am not sure whether it works if you choose 1
SWITCHES.CN.LinearCNratio=0;
% 
% 16. Turn on Denitrification and N20 flux from nitrification 
% Turn on == 1 or, Turn off =0
SWITCHES.CN.Denitrification=1;
% 
% 17. DEM_fraction_initialization
DEM_fraction_l_s_g_b_rate=[0 0 0 0];
Tot_UP_Amm_m2_print=zeros(nl_soil,1);
Tot_UP_Nit_m2_print=zeros(nl_soil,1);
CN_root=0;
% 
% 18. Urea initialization 
VARIABLES.Urea=0; % [g/m2]
SWITCHES.CN.Hydrolysis_on=1; % Hydrolysis switch
% 
% 19. Independent MIN and IMM : 0 Not use
SWITCHES.CN.Indep_MIN_IMM_on=0; 
% 
% 20. Poporato's fSd & fNd
SWITCHES.CN.PoporatofSdfNd_on=1; 
% 
% 20. Root exudation
SWITCHES.CN.exudation=1; 

% 21. Rainwater chemistry (10/22/2020 https://doi.org/10.1016/j.envres.2020.109872)
Rain_Nap = 3.05e-6; %mol/L = mol/kg
Rain_Ca2p = 4.7e-6; %mol/kg
Rain_Mg2p = 1.225e-6; %mol/kg
Rain_Kp = 7.3e-7; %mol/kg
Rain_Hp = 3.28978669598122E-06; %mol/kg TOTAL concentration 
Rain_Hpfert = 3.35689913429603E-05; %mol/kg TOTAL concentration for days when fert has been applied
Rain_Clm = 3.45e-6; %mol/kg
Rain_SO42m = 1.7e-5; %mol/kg
Rain_NO3m = 1.7e-5; %mol/kg
Rain_HCO3m = 8.5e-7; %mol/kg
Rain_O2 = 2.7e-4; %mol/kg
% 
% %-----------------------------------------------------------------------%%
% Initialization of Soil C & N                                            %
% %-----------------------------------------------------------------------%%
if from_previous_data ==0
%     1. Soil C, in litter, biomass, and humus pools [gC/m3]
%     Caution: used to calculate kl, kb, and kh, and also used to initially run the model.
%     Caution: kl, kb, and kh control the entire dinamics.
%     Caution: Not to put the values below 10
%     Cl and Ch is from EBI site.
% 
    VARIABLES.Cl =[2571.4;2564.1;2533.5;2375.05;2109.3;1628.5;977.35;579.8;357.9;158.65;....
        67.2;13.62];
    VARIABLES.Ch = [18523;18470;18249.8;17108.5;15194.5;11730.4;7040.5;4176.5;2578.5;...
        1143;484.1;98.2];
    VARIABLES.Cb = [358.595;357.575;353.3;331.225;294.16;227.1;136.305;80.86;49.915;...
        22.125;9.4;1.904];
%     
    figure()
        subplot(1,2,1)
            plot(rootfr,cumsum(dz))
            set(gca, 'YDir','reverse')
            ylabel('depth (m)')
            ylim([0 4.5])
            title('Root PDF')
        subplot(1,2,2)
            plot(VARIABLES.Cb,cumsum(dz))
            set(gca, 'YDir','reverse')
            ylabel('depth (m)')
            title('Initial Microbial Population')
%      close all
    if SWITCHES.CN_type == 0
        VARIABLES.Cl=VARIABLES.Cl+VARIABLES.Ch;
    end
%            
%     2. Soil N pool
%     Caution. Put any numbers. These will be adjusted after running
%     the first day [g/m3], VARIABLES.Nl and VARIABLES.CNl will be calculated
%     after average above-ground C/N ratio is calculated.
    VARIABLES.Nit=[2;2;2;2;2;2;2;2;2;2;2;2];
    Nit=VARIABLES.Nit;
%     Amm and NH3 will be recalculated in 'CN_compute_initial_Amm_NH3' with
%     pH value you set
    VARIABLES.Amm=[0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5];
    Amm=VARIABLES.Amm;
%     
    VARIABLES.Amm_mo=VARIABLES.Amm.*PARAMS.CN.a_Amm;
    Amm_mo=VARIABLES.Amm_mo;
    VARIABLES.Amm_im=VARIABLES.Amm.*(1-PARAMS.CN.a_Amm);
    Amm_im=VARIABLES.Amm_im;
    
    Amm_mo2im = zeros(nl_soil,1);
    Amm_im2mo = zeros(nl_soil,1);
    
    VARIABLES.NH3=VARIABLES.Amm.*0.05;
    NH3=VARIABLES.NH3;
%     
%  Initialize values in mol/kg for primary species for first day 
    VARIABLES.O2 = repmat(2.4813E-04,12,1); 
    O2_constant = repmat(2.4813E-04,12,1); % For now, holding this constant
    VARIABLES.CO2 = repmat(1.0606e-08,12,1); 
    VARIABLES.CO2gas = repmat(1.3456298E-06,12,1);
    VARIABLES.O2gas = repmat(0.116,12,1);
    VARIABLES.Hplus = repmat(3.28978669598122E-06,12,1); 
    VARIABLES.HCO3 = repmat(1.976E-14,12,1); 
    VARIABLES.CO32m = repmat(1.976E-14,12,1); 
    VARIABLES.N2 = repmat(1.0564e-08,12,1); 
    VARIABLES.glucose = repmat(0,12,1); 
    VARIABLES.cexg = repmat(0,12,1);
    VARIABLES.SOIL.cexgprev = zeros(12,365);
    VARIABLES.SOIL.cexfprev = zeros(12,365);
    VARIABLES.Cbaq = repmat(1E-12,12,1);
    VARIABLES.aceticacid = repmat(1,12,1);
    rainadd_Hp = zeros(12,365);
    rainadd_HCO3m = zeros(12,365);
    rainadd_NO3m = zeros(12,365);
    rainadd_Nap = zeros(12,365);
    rainadd_Ca2p = zeros(12,365);
    rainadd_Mg2p = zeros(12,365);
    rainadd_Kp = zeros(12,365);
    rainadd_O2 = zeros(12,365);
    VARIABLES.SiO2 = repmat(4.18006641001185E-04,12,1);
    VARIABLES.Al3p = repmat(5.52480585268130E-08,12,1);
    VARIABLES.Mg2p = repmat(1.96229806315877E-04,12,1);
    VARIABLES.Ca2p = repmat(1.38292166758517E-03,12,1);
    VARIABLES.Kp = repmat(2.22547613802969E-06,12,1);
    VARIABLES.Fe2p = repmat(7.55E-08,12,1);
    VARIABLES.Fe3p = repmat(6.39E-22,12,1);
    VARIABLES.TiO4H4 = repmat(2.76081894084137E-09,12,1);
    VARIABLES.Nap = repmat(7.12193976883976E-04,12,1);
%SOIL mineralogy
    VARIABLES.Albite = repmat(0.065384615,12,1); %bulk_surface_area 100
    VARIABLES.Anatase = repmat(0.003578947,12,1); %bulk_surface_area 100
    VARIABLES.Calcite = repmat(0.00903125,12,1); %bulk_surface_area 100
    VARIABLES.Chamosite = repmat(0.002833333,12,1); %specific_surface_area 1.0
    VARIABLES.Anthophyllite = repmat(0.006296296,12,1); %specific_surface_area 1.0
    VARIABLES.Kaolinite = repmat(0.006044444,12,1); %specific_surface_area 1.0
    VARIABLES.Microcline=  repmat(0.0756299212598425,12,1); %specific_surface_area 1.0
    VARIABLES.Muscovite = repmat(0.032785714,12,1); %specific_surface_area 1.0
    VARIABLES.Quartz = repmat(0.448415094,12,1); %bulk_surface_area 100
    VARIABLES.Goethite = repmat(0.010736842,12,1);% bulk_surface_area 100
    VARIABLES.Lepidocrocite= repmat(0.00799,12,1); % bulk_surface_area 100
    VARIABLES.Epidote = repmat(0.008742857,12,1); % bulk_surface_area 100
    VARIABLES.Magnetite = repmat(0.008286267,12,1); % bulk_surface_area 100
%     
end
% 
if Fix_Crop_Vaues==1
    adupsoil_layer_day_point=pre_adupsoil_layer_day_point;
end
% %-----------------------------------------------------------------------%%
%    PREALLOCATE VECTORS or MATRIXS
% %-----------------------------------------------------------------------%%          
CN_preallocation();
% 
% 
% %-----------------------------------------------------------------------%%
%    Load and organize parameters and forcings
% %-----------------------------------------------------------------------%%
CN_load_parameters_forcings(); 
% 
% 
% %-----------------------------------------------------------------------%%
%    Load forcings - for water age
% %-----------------------------------------------------------------------%%
CN_load_water_forcings()
% 
%     
% *************************************************************************%
% *************************************************************************%
% *                         Start the simulation                          *%
% *************************************************************************%            
% *************************************************************************%

%-----------------------------------------------------------------------%%
% Check How many year has been run
How_many_days_have_been_run=0;
Check_year=0;
for i=1:1:how_many_year
    detect_first_day_of_year(1)=1;
    detect_first_day_of_year(i+1)=1+365*i;
end
detect_first_day_of_year(end)=[];


%%-----------------------------------------------------------------------%%
disp('Simulation Start')

tic
% Start simulation - daily time step
for day=1:1:365*how_many_year
    if day <= 365
        doy(day) = day;
    else
        doy(day) = doy(day-365);
    end
     if  doy(day) == 239
         pause = 1;
     end
    
    
    %%-------------------------------------------------------------------%%
    % Check point to see how many days and years are being run
    for i=1:how_many_year
        if day == detect_first_day_of_year(i)
            Check_year=Check_year+1;
            % Check How many years has been run!
            Year=start_year+Check_year-1;
            Sim_Y_Dis=strcat('Year= ',num2str(Year));
            disp(Sim_Y_Dis)
        end
    end
    % Check How many days has been run!
    How_many_days_have_been_run=How_many_days_have_been_run+1;
        
 
    %%-------------------------------------------------------------------%%
    % Corn or soybean? Corn(1) or Soy(0) 
    if choose_crop == 2
        for j=1:size(corn_year,2)
            if Year == corn_year(j)
                iscorn=1;
                corn_days(day)=day;
            end
        end
        for j=1:size(soy_year,2)
            if Year == soy_year(j)
                iscorn=0;
                soy_days(day)=day;
            end
        end
    else
        iscorn=nan;
    end
    
    
    %%-------------------------------------------------------------------%%
    %susana edit
    % CR ratio and rootfr for corn & soybean
    dz = VERTSTRUC.dzsmm./1000; %depth of each layer in m
    depths = cumsum(dz);
    if choose_crop == 2
        if iscorn == 1
            if iscorn == 1
                te = 130; %(smith2013reduced + 3days)
                th = 302; %(smith2013reduced)
                VARIABLES.te = te;
                VARIABLES.th = th;
            end
            PARAMS.CN.CR_ratio = PARAMS.CR_ratio_Corn;
            HavestResidueRate = HavestResidueRate_Corn; 
            rootfr = [0.298643502; 0.182028575;0.093277394;0.060141897;...
                0.035497434;0.027550797;0.020486394;0.012263877;0.011239017;...
                0.006254449;0.003912412;0.003027905]; %based on method described in Amenu and Kumar;
        else
            te = 146; %(smith2013reduced + 3 days)
            th = 286;
            VARIABLES.te = te;
            VARIABLES.th = th;
            PARAMS.CN.CR_ratio = PARAMS.CR_ratio_Soy;
            HavestResidueRate = HavestResidueRate_Soy;
            rootfr = [0.064385184;0.336067495;0.285384533;0.148791167;...
                0.059117202;0.027340393;0.011278982;0.004275712;0.001988163;...
                0.000672498;0.000284008;0.000147718];
        end
    end

                
            
            
    
    %%-------------------------------------------------------------------%%
    
    % Put daily based varibales
    % - C/N ratio (for above ground)
    putCN=ALL_CNveg(day);
    for i=1:3
        CNveg(i)=putCN;
    end
    CNveg=CNveg(:);
    PARAMS.CN.CNveg=CNveg;
    
    % - C/N ratio (for above ground)
    putCN_root=ALL_CNveg_root(day);
    CNveg_root=putCN_root;
    PARAMS.CN.CNveg_root=CNveg_root;
    
    %%-------------------------------------------------------------------%%      
    % Hydrological variables
    % Caution: from PALMS, smp [MPa], layeruptake [mm/d], and qq [m/d]
    % Thay should be in [MPa], [m/d], and [m/d], repectively 
    % Because leaching hardly takes place when there is ice. It will be used in leaching.
    % sm: total water in soil (fluid water+soild water in soil). sm is used
    % when we compute the nitrogen in liquid or ice. 
    sm=vwc_layer_day_point(1:end-1,day)+ice_layer_day_point(1:end-1,day);
    smfull=vwc_layer_day_point(:,day)+ice_layer_day_point(:,day);
    real_sm=vwc_layer_day_point(:,day);
    
    ice=ice_layer_day_point(:,day);
    Ts=temp_layer_day_point(:,day);
    qq=qq_layer_day_point(:,day); 
    layeruptake_all=adupsoil_layer_day_point(:,day); %already in [m/d]
   
    % setting this as constant for now until 2-way coupled with MLCan
    VARIABLES.SOIL.hk = 14*10^-6; %m/s 
    
    % Rainfall Chemistry 
    if day> 1
        if PPT_mm(day) > 0
            sm_rain = vwc_layer_day_point(1:end-1,day) - vwc_layer_day_point(1:end-1,day-1); %approx diff in soil moisture due to rainfall
            sm_rain(sm_rain<0) = 0; %if negative, no rainwater added
            sm_rain_print(:,day) = sm_rain;
            if day > 125 && day < 130
                rainadd_Hp(:,day) = (sm_rain.*Rain_Hpfert)./vwc_layer_day_point(1:end-1,day); %mol/kg added 
            else
                rainadd_Hp(:,day) = (sm_rain.*Rain_Hp)./vwc_layer_day_point(1:end-1,day); %mol/kg added 
            end
            rainadd_HCO3m(:,day) = (sm_rain.*Rain_HCO3m)./vwc_layer_day_point(1:end-1,day); %mol/kg added 
            rainadd_NO3m(:,day) = (sm_rain.*Rain_NO3m)./vwc_layer_day_point(1:end-1,day); %mol/kg added 
            rainadd_Nap(:,day) = (sm_rain.*Rain_Nap)./vwc_layer_day_point(1:end-1,day); %mol/kg added 
            rainadd_Ca2p(:,day) = (sm_rain.*Rain_Ca2p)./vwc_layer_day_point(1:end-1,day); %mol/kg added 
            rainadd_Mg2p(:,day) = (sm_rain.*Rain_Mg2p)./vwc_layer_day_point(1:end-1,day); %mol/kg added 
            rainadd_Kp(:,day) = (sm_rain.*Rain_Kp)./vwc_layer_day_point(1:end-1,day); %mol/kg added 
            rainadd_O2(:,day) = (sm_rain.*Rain_O2)./vwc_layer_day_point(1:end-1,day); %mol/kg added 

        end
    end
    
   
    %%-------------------------------------------------------------------%%
    % Parameters and Variables    
    a_Amm = PARAMS.CN.a_Amm;%       a_Amm = fraction of dissolved ammonium [-]
    a_Nit = PARAMS.CN.a_Nit;%       a_Nit = fraction of dissolved nitrate [-]
    a_NH3 = PARAMS.CN.a_NH3;%       a_NH3 = fraction of dissolved ammonia [-]
    ki_Amm = PARAMS.CN.ki_Amm; % constant that determines the partition between NO3 and NH4 for immobilization
    ki_Nit = PARAMS.CN.ki_Nit; % constant that determines the partition between NO3 and NH4 for immobilization
    rr = PARAMS.CN.rr;%       rr = fraction of decomposed organic carbon that goes to respiration [-]
    Cl = VARIABLES.Cl;%       Cl = carbon concentration in litter pool [gC / m^3]
    Ch = VARIABLES.Ch;%       Ch = carbon concentration in humus pool [gC / m^3]
    Cb = VARIABLES.Cb;%       Cb = carbon concentration in biomass pool [gC / m^3]
    Nl = VARIABLES.Nl;%       Nl = nitrogen concentration in litter pool [gN / m^3]
    N2 = VARIABLES.N2; %        for Crunch
    CO2 = VARIABLES.CO2; %      for Crunch
    CO2gas = VARIABLES.CO2gas;          % for Crunch
    HCO3 = VARIABLES.HCO3; %        for Crunch
    CO32m = VARIABLES.CO32m; %      for Crunch
    O2 = VARIABLES.O2; %            for Crunch
    O2gas = VARIABLES.O2gas; %              for Crunch
    Hplus = VARIABLES.Hplus; %          for Crunch 
    pH = VARIABLES.pH;
    glucose = VARIABLES.glucose; % concentration of glucose left in soil g/m3
    Cbaq =  VARIABLES.Cbaq;
    aceticacid =  VARIABLES.aceticacid;
    aceticacidoverride = zeros(size(Cb));
    cnmodeldCbdt = zeros(12,365);
    cbpassedtocrunch = zeros(12,365);
    SiO2 =  VARIABLES.SiO2;
    Al3p =  VARIABLES.Al3p;
    Mg2p =  VARIABLES.Mg2p;
    Ca2p =  VARIABLES.Ca2p;
    Kp =  VARIABLES.Kp;
    Fe2p =  VARIABLES.Fe2p;
    Fe3p =  VARIABLES.Fe3p;
    TiO4H4 =  VARIABLES.TiO4H4;
    Nap =  VARIABLES.Nap;
    cexg = VARIABLES.cexg;
    Albite = VARIABLES.Albite;
    Anatase = VARIABLES.Anatase;
    Calcite = VARIABLES.Calcite;
    Chamosite = VARIABLES.Chamosite;
    Anthophyllite = VARIABLES.Anthophyllite;
    Kaolinite = VARIABLES.Kaolinite;
    Microcline = VARIABLES.Microcline;
    Muscovite = VARIABLES.Muscovite;
    Quartz = VARIABLES.Quartz;
    Goethite = VARIABLES.Goethite;
    Lepidocrocite = VARIABLES.Lepidocrocite;
    Epidote = VARIABLES.Epidote;
    Magnetite = VARIABLES.Magnetite;
                
	if length(Nl)> 12
        Nl = Nl(2:end);
	end
    Amm = VARIABLES.Amm;%       Amm = ammonium concentration in soil [gN / m^3]

    %%-------------------------------
    % NEED CHANGING
    Amm_mo = VARIABLES.Amm_mo;%       Mobile Amm [gN / m^3]
    Amm_im = VARIABLES.Amm_im;%       Immobile Amm [gN / m^3]
    
    %%-------------------------------
    
    Nit = VARIABLES.Nit;%       Nit = nitrate concentration in soil [gN / m^3]
    NH3 = VARIABLES.NH3;      % NH3 = Ammonia concentration in soil [gN/m3]
 
    
    dtime = CONSTANTS.dtime;%       dtime = Time step [1800 s]
    nl_soil=PARAMS.nl_soil;%       nl_soil = # soil layers
    dz_mm = VERTSTRUC.dzsmm;%       dz_mm = grid depth [mm]
    nspecies = PARAMS.CanStruc.nspecies; % number of species
     

    
%**************************************************************************
%    PREALLOCATE VECTORS 
%**************************************************************************           
    PHI = zeros(nl_soil,1);
    phi = ones(nl_soil,1);  


%**************************************************************************
%    COMPUTE ENVIRONMENTAL FACTORS 
%**************************************************************************            
    CN_envfactor (); 
    
   
%*************************************************************************
%      Root exudation
%************ *************************************************************     
    if day == 1
        Brprev = 0;
        VARIABLES.SOIL.Br = Brprev;
    end


    if SWITCHES.CN.exudation
        [gmu, Cbglu, Gtot,Ftot, fflav,cexf,cexg, Bp, Br,VARIABLES] = ...
            CN_exudate(doy,day,iscorn,VARIABLES,VERTSTRUC,qq,layeruptake_all,...
            PARAMS,smfull,rootfr,BulkDensity,te,th);
        VARIABLES.SOIL.Gtot = Gtot;
        VARIABLES.SOIL.Ftot = Ftot;
        VARIABLES.SOIL.fflav = fflav;
        VARIABLES.SOIL.Br = Br;          
    else 
        Cbglu = [];
        cexg
    end
    
%% CRUNCH/TRANSFORM FIRST
%*************************************************************************
%      Crunch
%************ *************************************************************   
% Denitrification rate factors: TD.*WD.*MD *rate
    [F_D] = CN_denitrification_factor(VARIABLES, VERTSTRUC, smp, temp_layer_day_point(:,day), Cb);
% Nitrification rate factor: (1-fflav)*fSn.*fTn * rate
    [F_N] = CN_nitrification_factor(fSn,fTn,fflav);

%declare var name as they appear in crunch
vartocrunch = {'H+';'C5H7O2NO2(aq)';'Acetic_acid(aq)';'C6H12O6';...
    'NH4+';'NO3-';'N2(aq)';'O2(aq)';'CO2(aq)';'SiO2(aq)';'Al+++';...
    'Mg++';'Ca++';'K+';'Fe++';'Fe+++';'Ti(OH)4(aq)';'Na+'};   
mincrunch = {'C5H7O2NO2(s)';'Albite';'Anatase';'Calcite';'Chamosite';...
    'Anthophyllite';'Kaolinite';'Microcline';'Muscovite';'Quartz';...
    'Goethite';'Lepidocrocite';'Epidote';'Magnetite'};   

for ss=1:N
    if Ts(ss) > 0
        %prepare string overwrite for primary species + concentration
        primarytocrunch(1) = strcat(vartocrunch{1},{' '},num2str(Hplus(ss))); %mol/kg
        primarytocrunch(2) = strcat(vartocrunch{2},{' '},num2str(Cbaq(ss))); % mol/kg
        primarytocrunch(3) = strcat(vartocrunch{3},{' '},num2str(aceticacid(ss))); % mol/kg

        totglucose = (cexg(ss) + glucose(ss))./180.156.*0.001./sm(ss);
        primarytocrunch(4) = strcat(vartocrunch{4},{' '},num2str(totglucose)); %g/m3 converted to mol/kg

        totNH4 = (NH3(ss))./17.03052.*0.001./sm(ss) +...
                (Amm(ss))./18.03846.*0.001./sm(ss); 
            if totNH4 < 10^-12
            nh4crunchin = 10^-12;
        else
            nh4crunchin = totNH4 ;
        end
        primarytocrunch(5) = strcat(vartocrunch{5},{' '},num2str(nh4crunchin)); %mol/kg

        totNit = (Nit(ss))./62.0049.*0.001./sm(ss);
        primarytocrunch(6) = strcat(vartocrunch{6},{' '},num2str(totNit)); %g/m3 converted to mol/kg

        primarytocrunch(7) = strcat(vartocrunch{7},{' '},num2str(N2(ss))); %mol/kg

        primarytocrunch(8) = strcat(vartocrunch{8},{' '},num2str(O2_constant(ss))); %mol/kg

        totCO2 = CO2(ss) + HCO3(ss) + CO32m(ss);
        primarytocrunch(9) = strcat(vartocrunch{9},{' '},num2str(totCO2)); %mol/kg

        primarytocrunch(10) = strcat(vartocrunch{10},{' '},num2str(SiO2(ss))); %mol/kg
        primarytocrunch(11) = strcat(vartocrunch{11},{' '},num2str(Al3p(ss))); %mol/kg
        primarytocrunch(12) = strcat(vartocrunch{12},{' '},num2str(Mg2p(ss))); %mol/kg
        primarytocrunch(13) = strcat(vartocrunch{13},{' '},num2str(Ca2p(ss))); %mol/kg
        primarytocrunch(14) = strcat(vartocrunch{14},{' '},num2str(Kp(ss))); %mol/kg
        primarytocrunch(15) = strcat(vartocrunch{15},{' '},num2str(Fe2p(ss))); %mol/kg
        primarytocrunch(16) = strcat(vartocrunch{16},{' '},num2str(Fe3p(ss))); %mol/kg
        primarytocrunch(17) = strcat(vartocrunch{17},{' '},num2str(TiO4H4(ss))); %mol/kg
        primarytocrunch(18) = strcat(vartocrunch{18},{' '},num2str(Nap(ss))); %mol/kg

        convertcb = (1/12.0107)*(1/5)*(1/113)*(3125)*(100^(-3));
        mintocrunch(1) = strcat(mincrunch(1),{' '},num2str((Cb(ss)).*convertcb)); %m3/m3         
        mintocrunch(2) = strcat(mincrunch(2),{' '},num2str(Albite(ss)),{' '},'bulk_surface_area',{' '}, num2str(100)); %m3/m3         
        mintocrunch(3) = strcat(mincrunch(3),{' '},num2str(Anatase(ss)),{' '},'bulk_surface_area',{' '}, num2str(100)); %m3/m3         
        mintocrunch(4) = strcat(mincrunch(4),{' '},num2str(Calcite(ss)),{' '},'bulk_surface_area',{' '}, num2str(100)); %m3/m3         
        mintocrunch(5) = strcat(mincrunch(5),{' '},num2str(Chamosite(ss)),{' '},'specific_surface_area',{' '}, num2str(1.0)); %m3/m3         
        mintocrunch(6) = strcat(mincrunch(6),{' '},num2str(Anthophyllite(ss)),{' '},'specific_surface_area',{' '}, num2str(1.0)); %m3/m3         
        mintocrunch(7) = strcat(mincrunch(7),{' '},num2str(Kaolinite(ss)),{' '},'specific_surface_area',{' '}, num2str(1.0)); %m3/m3         
        mintocrunch(8) = strcat(mincrunch(8),{' '},num2str(Microcline(ss)),{' '},'specific_surface_area',{' '}, num2str(1.0)); %m3/m3         
        mintocrunch(9) = strcat(mincrunch(9),{' '},num2str(Muscovite(ss)),{' '},'specific_surface_area',{' '}, num2str(1.0)); %m3/m3         
        mintocrunch(10) = strcat(mincrunch(10),{' '},num2str(Quartz(ss)),{' '},'bulk_surface_area',{' '}, num2str(100)); %m3/m3         
        mintocrunch(11) = strcat(mincrunch(11),{' '},num2str(Goethite(ss)),{' '},'bulk_surface_area',{' '}, num2str(100)); %m3/m3         
        mintocrunch(12) = strcat(mincrunch(12),{' '},num2str(Lepidocrocite(ss)),{' '},'bulk_surface_area',{' '}, num2str(100)); %m3/m3         
        mintocrunch(13) = strcat(mincrunch(13),{' '},num2str(Epidote(ss)),{' '},'bulk_surface_area',{' '}, num2str(100)); %m3/m3         
        mintocrunch(14) = strcat(mincrunch(14),{' '},num2str(Magnetite(ss)),{' '},'bulk_surface_area',{' '}, num2str(100)); %m3/m3         

        %Get soil temp
        soiltemp = strcat('temperature',{' '},num2str(temp_layer_day_point(ss,day))); %deg C

        %Prepare output files
        outtocrunch(primarytocrunch,mintocrunch,soiltemp,F_N(ss),F_D(ss),dz(ss),sm(ss));
        fclose('all');
        %Run crunch
        system('CrunchTope.exe');
        fclose('all');
        %Get state updates
        [pH_crunch,Cb_aq,aceticacid_crunch,C6H12O6_crunch,...
                    NH4_crunch,NO3_crunch,N2_crunch,O2_crunch,...
                    CO2_crunch,HCO3_crunch,CO32m_crunch,Hplus_crunch,...
                    Cb_crunch,NH3_crunch,CO2gas_crunch,O2gas_crunch,...
                    SiO2_crunch,Al3p_crunch,Mg2p_crunch,Ca2p_crunch,...
                    Kp_crunch,Fe2p_crunch,Fe3p_crunch,TiO4H4_crunch,Nap_crunch,...
                    Albite_crunch, Anatase_crunch, Calcite_crunch,Chamosite_crunch,...
                    Anthophyllite_crunch, Kaolinite_crunch, Microcline_crunch,...
                    Muscovite_crunch, Quartz_crunch,Goethite_crunch,...
                    Lepidocrocite_crunch,Epidote_crunch, Magnetite_crunch] = backfromcrunch();
        pH_cr(ss) = pH_crunch;
        Hplus_cr(ss) = Hplus_crunch;
        Cb_cr(ss) = Cb_crunch./convertcb;
        Cbaq_cr(ss) = Cb_aq; 
        glucose_cr(ss) = C6H12O6_crunch.*180.156./0.001.*sm(ss); %back to g/m3
        NH3_cr(ss) = NH3_crunch.*17.03052./0.001.*sm(ss); 
        Amm_cr(ss) = NH4_crunch.*18.03846./0.001.*sm(ss); 
        N2_cr(ss) = N2_crunch;
        Nit_cr(ss) = NO3_crunch.*62.0049./0.001.*sm(ss);
        O2_cr(ss) = O2_crunch;
        O2gas_cr(ss) = O2gas_crunch;
        CO2_cr(ss) = CO2_crunch;   
        CO2gas_cr(ss) = CO2gas_crunch;
        HCO3_cr(ss) = HCO3_crunch;   
        CO32m_cr(ss) = CO32m_crunch;
        aceticacid_cr(ss) = aceticacid_crunch;
        SiO2_cr(ss) = SiO2_crunch;
        Al3p_cr(ss) = Al3p_crunch;
        Mg2p_cr(ss) = Mg2p_crunch;
        Ca2p_cr(ss) = Ca2p_crunch;
        Kp_cr(ss) = Kp_crunch;
        Fe2p_cr(ss) = Fe2p_crunch;
        Fe3p_cr(ss) = Fe3p_crunch;
        TiO4H4_cr(ss) = TiO4H4_crunch;
        Nap_cr(ss) = Nap_crunch;
        Albite_cr(ss) = Albite_crunch;
        Anatase_cr(ss) = Anatase_crunch;
        Calcite_cr(ss) = Calcite_crunch;
        Chamosite_cr(ss) = Chamosite_crunch;
        Anthophyllite_cr(ss) = Anthophyllite_crunch;
        Kaolinite_cr(ss) = Kaolinite_crunch;
        Microcline_cr(ss) = Microcline_crunch;
        Muscovite_cr(ss) = Muscovite_crunch;
        Quartz_cr(ss) = Quartz_crunch;
        Goethite_cr(ss) = Goethite_crunch;
        Lepidocrocite_cr(ss) = Lepidocrocite_crunch;
        Epidote_cr(ss) = Epidote_crunch;
        Magnetite_cr(ss) = Magnetite_crunch;
        close all;
    else
        pH_cr(ss) = pH(ss);
        Hplus_cr(ss) = Hplus(ss);
        Cb_cr(ss) = Cb(ss);
        Cbaq_cr(ss) = Cbaq(ss); 
        glucose_cr(ss) = cexg(ss) + glucose(ss); %back to g/m3
        NH3_cr(ss) = NH3(ss); 
        Amm_cr(ss) = Amm(ss); 
        N2_cr(ss) = N2(ss);
        Nit_cr(ss) = Nit(ss);
        O2_cr(ss) = O2(ss);
        O2gas_cr(ss) = O2gas(ss);
        CO2_cr(ss) = CO2(ss);
        CO2gas_cr(ss) = CO2gas(ss);
        HCO3_cr(ss) = HCO3(ss);   
        CO32m_cr(ss) = CO32m(ss);
        aceticacid_cr(ss) = aceticacid(ss);
        SiO2_cr(ss) = SiO2(ss);
        Al3p_cr(ss) = Al3p(ss);
        Mg2p_cr(ss) = Mg2p(ss);
        Ca2p_cr(ss) = Ca2p(ss);
        Kp_cr(ss) = Kp(ss);
        Fe2p_cr(ss) = Fe2p(ss);
        Fe3p_cr(ss) = Fe3p(ss);
        TiO4H4_cr(ss) = TiO4H4(ss);
        Nap_cr(ss) = Nap(ss);
        Albite_cr(ss) = Albite(ss);
        Anatase_cr(ss) = Anatase(ss);
        Calcite_cr(ss) = Calcite(ss);
        Chamosite_cr(ss) = Chamosite(ss);
        Anthophyllite_cr(ss) = Anthophyllite(ss);
        Kaolinite_cr(ss) = Kaolinite(ss);
        Microcline_cr(ss) = Microcline(ss);
        Muscovite_cr(ss) = Muscovite(ss);
        Quartz_cr(ss) = Quartz(ss);
        Goethite_cr(ss) = Goethite(ss);
        Lepidocrocite_cr(ss) = Lepidocrocite(ss);
        Epidote_cr(ss) = Epidote(ss);
        Magnetite_cr(ss) = Magnetite(ss);
    end
end

pH_cr = reshapearray(pH_cr);
Hplus_cr = reshapearray(Hplus_cr);
Cb_cr = reshapearray(Cb_cr);
glucose_cr = reshapearray(glucose_cr);
NH3_cr = reshapearray(NH3_cr);
Amm_cr = reshapearray(Amm_cr);
Nit_cr = reshapearray(Nit_cr);
O2_cr = reshapearray(O2_cr);
O2gas_cr = reshapearray(O2gas_cr);
CO2_cr = reshapearray(CO2_cr);
CO2gas_cr = reshapearray(CO2gas_cr);
N2_cr = reshapearray(N2_cr);
HCO3_cr = reshapearray(HCO3_cr);
CO32m_cr = reshapearray(CO32m_cr);
Cbaq_cr = reshapearray(Cbaq_cr);
aceticacid_cr = reshapearray(aceticacid_cr);
SiO2_cr = reshapearray(SiO2_cr);
Al3p_cr = reshapearray(Al3p_cr);
Mg2p_cr = reshapearray(Mg2p_cr);
Ca2p_cr = reshapearray(Ca2p_cr);
Kp_cr = reshapearray(Kp_cr);
Fe2p_cr = reshapearray(Fe2p_cr);
Fe3p_cr = reshapearray(Fe3p_cr);
TiO4H4_cr = reshapearray(TiO4H4_cr);
Nap_cr = reshapearray(Nap_cr);
Albite_cr = reshapearray(Albite_cr);
Calcite_cr = reshapearray(Calcite_cr);
Chamosite_cr = reshapearray(Chamosite_cr);
Anthophyllite_cr = reshapearray(Anthophyllite_cr);
Kaolinite_cr = reshapearray(Kaolinite_cr);
Microcline_cr = reshapearray(Microcline_cr);
Muscovite_cr = reshapearray(Muscovite_cr);
Quartz_cr = reshapearray(Quartz_cr);
Goethite_cr = reshapearray(Goethite_cr);
Lepidocrocite_cr = reshapearray(Lepidocrocite_cr);
Epidote_cr = reshapearray(Epidote_cr);
Magnetite_cr = reshapearray(Magnetite_cr);

%**************************************************************************
%    BIOTURBATION 
%**************************************************************************
    if SWITCHES.CN.Bioturbation
        if SWITCHES.CN.Bioturbation_new == 0
            [Cl, VARIABLES] = CN_bioturbation (PARAMS, VARIABLES, CONSTANTS,...
                FORCING, VERTSTRUC, Cl, day);
        elseif SWITCHES.CN.Bioturbation_new == 1
            [VARIABLES] = CN_bioturbation_new (PARAMS, VARIABLES, CONSTANTS,...
                FORCING, VERTSTRUC, SWITCHES, Cl, Cb, Ch, day);            
        end
    end  
     
   
%*************************************************************************
%     COMPUTE RATES OF CHANGES IN POOLS CONCENTRATION. (LITTER, HUMUS, BIOMASS) 
%*************************************************************************  
    [CNa, ADD, CNo, OUT, ADD_bio, ADD_ex, ADD_net, Before_Bioturbation_ADD] = CN_addlitter (FORCING, PARAMS,...
        VERTSTRUC, VARIABLES, CONSTANTS, SWITCHES, rootfr, choose_crop,...
        TBMLaNoLitterBack, RootDeathWhenHarvest, HavestResidueRate, iscorn, day);

    
%**************************************************************************
%        Atmospheric deposition
%************************************************************************** 

    Depo_Nit_m3=zeros (nl_soil,1);
    Depo_Amm_m3=zeros (nl_soil,1);
    if SWITCHES.CN.Ndeposition == 1
        Depo_Nit_m3(1)=Deposition_NO3_N_Apply(day)./(dz_mm(1)./1000);
        Depo_Amm_m3(1)=Deposition_NH4_N_Apply(day)./(dz_mm(1)./1000);       
    end
        Depo_Nit_m2=Depo_Nit_m3.*(dz_mm./1000);
        Depo_Amm_m2=Depo_Amm_m3.*(dz_mm./1000);

    
%**************************************************************************
%        Fertilizer application
%**************************************************************************  
    if CORN == 1
        [Fert_Amm_m3, Fert_Nit_m3, Ins_Vol_m2_1st_layer, NH3_Urea_m2, NH4_Urea_m2, VARIABLES, SWITCHES]...
            = CN_fertilizer(VARIABLES, VERTSTRUC, SWITCHES, day, IN_Nit_Fertilizer,...
            IN_Amm_Fertilizer, IN_Organic_Fertilizer, nl_soil, PPT_mm, Wind_ms, BulkDensity,...
            sm, Ts, smp, WFPS, Clay, pH);
    %       
        Ins_Vol_m2=zeros(nl_soil,1);
        Ins_Vol_m2(1)=Ins_Vol_m2_1st_layer;
        Ins_Vol_m3=Ins_Vol_m2./(dz_mm./1000);
        NH3_Urea_m3=NH3_Urea_m2./(dz_mm./1000);
        NH4_Urea_m3=NH4_Urea_m2./(dz_mm./1000);

        Fert_Nit_m2=Fert_Nit_m3.*(dz_mm./1000);
        Fert_Amm_m2=Fert_Amm_m3.*(dz_mm./1000);

        Urea(1,1) = VARIABLES.Urea;
        Urea(2:length(Cl),1) = 0;
    else
        Ins_Vol_m2=zeros(nl_soil,1);
        NH3_Urea_m2=zeros(nl_soil,1);
        NH4_Urea_m2=zeros(nl_soil,1);
        Fert_Nit_m3=zeros(nl_soil,1);
        Fert_Amm_m3=zeros(nl_soil,1);

        Ins_Vol_m3=Ins_Vol_m2./(dz_mm./1000);
        NH3_Urea_m3=NH3_Urea_m2./(dz_mm./1000).*0;
        NH4_Urea_m3=NH4_Urea_m2./(dz_mm./1000).*0;

        Fert_Nit_m2=Fert_Nit_m3.*(dz_mm./1000).*0;
        Fert_Amm_m2=Fert_Amm_m3.*(dz_mm./1000).*0;

        Urea(1,1) = VARIABLES.Urea.*0;
        Urea(2:length(Cl),1) = 0;
    end
     
% %**************************************************************************
% %       Volatilization
% %**************************************************************************
    % Volatlization is gas exchange between soil surface and atmosphere
    Volat= zeros(nl_soil,1);
    NH3_m2=NH3_cr(1).*(dz_mm(1)/1000); %susana - NH3_cr not NH3
    
    [Volat_m2] = CN_volatilization(NH3_m2, Ts, WFPS, Clay);
    Volat(1)=Volat_m2/(dz_mm(1)/1000);
% %**************************************************************************
% %       Aqueous/gaseous CO2 loss
% %**************************************************************************
   %Aqueous CO2 building up; should be lost as gas
   alpha = 0.1; %fraction of CO2 that remains aqueous
   resp2aq = alpha.*CO2_cr;
   resp2gas = (1-alpha).*CO2_cr;
        
%**************************************************************************
%     COMPUTE NITRATE AND AMMONIUM WATER UPTAKE FLUX [gN / m^2 / d]
%**************************************************************************            
[SUP_amm_m2, SUP_nit_m2, SUP_amm_all_m2, SUP_nit_all_m2, SUP_amm_active_m2, SUP_nit_active_m2, SUP_amm_passive_m2, SUP_nit_passive_m2,...
        Amm_DEM, Nit_DEM, Fix_Nit_m2, Fix_Amm_m2, DEM_fraction_l_s_g_b, VARIABLES] = ...
        CN_nuptake (VARIABLES, PARAMS, SWITCHES, VERTSTRUC,...  
        a_Nit, a_Amm, sm, layeruptake_all, CARBONS, CNratio_leaf, CNratio_stem, CNratio_grain, CNratio_root, rootfr, choose_crop,...
        scenarios, totlail_layer_day_point, Add_Carbon, how_many_year,...
        detect_first_day_of_year, day, iscorn,...
        DEM_fraction_l_s_g_b_rate,Check_year,Tot_UP_Amm_m2_print,Tot_UP_Nit_m2_print,Amm_cr,Nit_cr);

    
%     Compute uptake in units of [gr / m^3 / d]
%     N uptake from soil
    SUP_nit = SUP_nit_m2./(dz_mm/1000);   % [gr/m^3/d]
    SUP_amm = SUP_amm_m2./(dz_mm/1000);   % [gr/m^3/d]
        
    SUP_nit_all = SUP_nit_all_m2./repmat((dz_mm/1000),1,nspecies); % [gr/m^3/d]        
    SUP_amm_all = SUP_amm_all_m2./repmat((dz_mm/1000),1,nspecies); % [gr/m^3/d]

    SUP_nit_active=SUP_nit_active_m2./(dz_mm/1000);
    SUP_amm_active=SUP_amm_active_m2./(dz_mm/1000);
    
    SUP_nit_passive=SUP_nit_passive_m2./(dz_mm/1000);
    SUP_amm_passive=SUP_amm_passive_m2./(dz_mm/1000);
    
%     N uptake from N fixation
    Fix_Nit=Fix_Nit_m2./(dz_mm/1000);   % [gr/m^3/d]
    Fix_Amm=Fix_Amm_m2./(dz_mm/1000);   % [gr/m^3/d]
    
%     N uptake from N remobilization
    Nre_Nit = zeros(size(Fix_Nit));
    Nre_Amm = zeros(size(Fix_Amm));
    
%     N uptake from soil, fixation, and n remobilization
    Tot_UP_Nit=SUP_nit+Fix_Nit+Nre_Nit; % [gr/m^3/d]
    Tot_UP_Amm=SUP_amm+Fix_Amm+Nre_Amm; % [gr/m^3/d]
    
%     Compute total Nitrogen Uptake [gr/m2/d].        
    Tot_UP_N_m2 = sum(Tot_UP_Amm.*(dz_mm/1000)) + sum(Tot_UP_Nit.*(dz_mm/1000));   %[gr/m^2/d]
    Tot_SUP_N_m2 = sum(SUP_amm.*(dz_mm/1000)) + sum(SUP_nit.*(dz_mm/1000));   %[gr/m^2/d]
    
%     Check the mass balance
    mberracNUP=Tot_UP_N_m2-sum(SUP_nit_active_m2)-sum(SUP_nit_passive_m2)-...
        sum(SUP_amm_active_m2)-sum(SUP_amm_passive_m2)-sum(Fix_Nit_m2)-sum(Fix_Amm_m2);

% %**************************************************************************
% %     COMPUTE RATES OF CHANGES IN POOLS CONCENTRATION. (LITTER, HUMUS, BIOMASS) 
% %**************************************************************************  
    [dCl_dt, dNl_dt, dCh_dt, dNh_dt, dCb_dt, dNb_dt, DECl, VARIABLES, DECh, Nladd, rh] =...
        CN_dynamics_new (VARIABLES, PARAMS, SWITCHES, fSd, fTd, phi, PHI, ADD, CNa, OUT, CNo);    
    respC = (DECl+DECh).*dz*rr;   %[grC/m2/d] Respiration
    respCgas = (1-alpha).*respC;
    respCaqueous = alpha.*respC;
    
%% TRANSPORT   
%**************************************************************************
%     COMPUTE NITRATE AND AMMONIUM LEACHING FLUX [gN / m^3 / d]
%************************************************************************** 
% UPDATE CONCENTRATIONS FIRST, THEN LEACH
omega = 10^-13;
%     NH3(aq) 
        dNH3_dt_gl = NH3_Urea_m3 - Volat ;
        NH3_gl = NH3_cr + dNH3_dt_gl;
        NH3_gl(abs(NH3_gl)<0)=omega; 
               
%     NH4+
        dAmm_dt_gl =  Depo_Amm_m3 + Fert_Amm_m3 + NH4_Urea_m3 - SUP_amm;
        Amm_gl = Amm_cr + dAmm_dt_gl;
        % To avoid non-physical negative values
        Amm_gl(Amm_gl<0) = 1e-10;

%     NO3-
        dNit_dt_gl =  Depo_Nit_m3 + Fert_Nit_m3 - SUP_nit;% ;
        Nit_gl = Nit_cr + dNit_dt_gl;
        % To avoid non-physical negative values
        Nit_gl(Nit_gl<0) = 1e-10;

%     CO2
        dCO2_dt_gl = -resp2gas + respCaqueous; %respCaqueous;
        CO2_gl = CO2_cr + dCO2_dt_gl;
        CO2_gl(abs(CO2_gl)<0)=0;  

% Primary species from Crunch that can be transported must be
% transported here. 
    [LCH_Hplus_m2, TLCH_Hplus_m2, flow_Hplus_m2, LCH_Hplus_Top_Up_m2, LCH_Hplus_Bottom_Down_m2, LCH_Hplus_Top_Down_m2, LCH_Hplus_Bottom_Up_m2]...
         = CN_leach(PARAMS, Hplus_cr, 1, qq, smfull, dz_mm);   
    [LCH_glucose_m2, TLCH_glucose_m2, flow_glucose_m2, LCH_glucose_Top_Up_m2, LCH_glucose_Bottom_Down_m2, LCH_glucose_Top_Down_m2, LCH_glucose_Bottom_Up_m2]...
        = CN_leach(PARAMS, glucose_cr, 1, qq, smfull, dz_mm);  
    [LCH_NH3_m2, TLCH_NH3_m2, flow_NH3_m2, LCH_NH3_Top_Up_m2, LCH_NH3_Bottom_Down_m2, LCH_NH3_Top_Down_m2, LCH_NH3_Bottom_Up_m2]...
        = CN_leach(PARAMS, NH3_gl, 1, qq, smfull, dz_mm);
    [LCH_Amm_m2, TLCH_Amm_m2, flow_Amm_m2, LCH_Amm_Top_Up_m2, LCH_Amm_Bottom_Down_m2, LCH_Amm_Top_Down_m2, LCH_Amm_Bottom_Up_m2]...
        = CN_leach(PARAMS, Amm_gl.*a_Amm, 1, qq, smfull, dz_mm);
    [LCH_Nit_m2, TLCH_Nit_m2, flow_Nit_m2, LCH_Nit_Top_Up_m2, LCH_Nit_Bottom_Down_m2, LCH_Nit_Top_Down_m2, LCH_Nit_Bottom_Up_m2]...
        = CN_leach(PARAMS, Nit_gl, 1, qq, smfull, dz_mm);   
    [LCH_O2_m2, TLCH_O2_m2, flow_O2_m2, LCH_O2_Top_Up_m2, LCH_O2_Bottom_Down_m2, LCH_O2_Top_Down_m2, LCH_O2_Bottom_Up_m2]...
        = CN_leach(PARAMS, O2_cr, 1, qq, smfull, dz_mm); 
    [LCH_CO2_m2, TLCH_CO2_m2, flow_CO2_m2, LCH_CO2_Top_Up_m2, LCH_CO2_Bottom_Down_m2, LCH_CO2_Top_Down_m2, LCH_CO2_Bottom_Up_m2]...
        = CN_leach(PARAMS, CO2_gl, 1, qq, smfull, dz_mm); 
    [LCH_HCO3_m2, TLCH_HCO3_m2, flow_HCO3_m2, LCH_HCO3_Top_Up_m2, LCH_HCO3_Bottom_Down_m2, LCH_HCO3_Top_Down_m2, LCH_HCO3_Bottom_Up_m2]...
        = CN_leach(PARAMS, HCO3_cr, 1, qq, smfull, dz_mm);
    [LCH_CO32m_m2, TLCH_CO32m_m2, flow_CO32m_m2, LCH_CO32m_Top_Up_m2, LCH_CO32m_Bottom_Down_m2, LCH_CO32m_Top_Down_m2, LCH_CO32m_Bottom_Up_m2]...
        = CN_leach(PARAMS, CO32m_cr, 1, qq, smfull, dz_mm);
    [LCH_N2_m2, TLCH_N2_m2, flow_N2_m2, LCH_N2_Top_Up_m2, LCH_N2_Bottom_Down_m2, LCH_N2_Top_Down_m2, LCH_N2_Bottom_Up_m2]...
        = CN_leach(PARAMS, N2_cr, 1, qq, smfull, dz_mm);
    [LCH_Cbaq_m2, TLCH_Cbaq_m2, flow_Cbaq_m2, LCH_Cbaq_Top_Up_m2, LCH_Cbaq_Bottom_Down_m2, LCH_Cbaq_Top_Down_m2, LCH_Cbaq_Bottom_Up_m2]...
        = CN_leach(PARAMS, Cbaq_cr, 1, qq, smfull, dz_mm);
    [LCH_aceticacid_m2, TLCH_aceticacid_m2, flow_aceticacid_m2, LCH_aceticacid_Top_Up_m2, LCH_aceticacid_Bottom_Down_m2, LCH_aceticacid_Top_Down_m2, LCH_aceticacid_Bottom_Up_m2]...
        = CN_leach(PARAMS, aceticacid_cr, 1, qq, smfull, dz_mm);
    [LCH_SiO2_m2, TLCH_SiO2_m2, flow_SiO2_m2, LCH_SiO2_Top_Up_m2, LCH_SiO2_Bottom_Down_m2, LCH_SiO2_Top_Down_m2, LCH_SiO2_Bottom_Up_m2]...
        = CN_leach(PARAMS, SiO2_cr, 1, qq, smfull, dz_mm);
   [LCH_Al3p_m2, TLCH_Al3p_m2, flow_Al3p_m2, LCH_Al3p_Top_Up_m2, LCH_Al3p_Bottom_Down_m2, LCH_Al3p_Top_Down_m2, LCH_Al3p_Bottom_Up_m2]...
        = CN_leach(PARAMS, Al3p_cr, 1, qq, smfull, dz_mm);
   [LCH_Mg2p_m2, TLCH_Mg2p_m2, flow_Mg2p_m2, LCH_Mg2p_Top_Up_m2, LCH_Mg2p_Bottom_Down_m2, LCH_Mg2p_Top_Down_m2, LCH_Mg2p_Bottom_Up_m2]...
        = CN_leach(PARAMS, Mg2p_cr, 1, qq, smfull, dz_mm);
    [LCH_Ca2p_m2, TLCH_Ca2p_m2, flow_Ca2p_m2, LCH_Ca2p_Top_Up_m2, LCH_Ca2p_Bottom_Down_m2, LCH_Ca2p_Top_Down_m2, LCH_Ca2p_Bottom_Up_m2]...
        = CN_leach(PARAMS, Ca2p_cr, 1, qq, smfull, dz_mm);
    [LCH_Kp_m2, TLCH_Kp_m2, flow_Kp_m2, LCH_Kp_Top_Up_m2, LCH_Kp_Bottom_Down_m2, LCH_Kp_Top_Down_m2, LCH_Kp_Bottom_Up_m2]...
        = CN_leach(PARAMS, Kp_cr, 1, qq, smfull, dz_mm);
    [LCH_Fe2p_m2, TLCH_Fe2p_m2, flow_Fe2p_m2, LCH_Fe2p_Top_Up_m2, LCH_Fe2p_Bottom_Down_m2, LCH_Fe2p_Top_Down_m2, LCH_Fe2p_Bottom_Up_m2]...
        = CN_leach(PARAMS, Fe2p_cr, 1, qq, smfull, dz_mm);
    [LCH_Fe3p_m2, TLCH_Fe3p_m2, flow_Fe3p_m2, LCH_Fe3p_Top_Up_m2, LCH_Fe3p_Bottom_Down_m2, LCH_Fe3p_Top_Down_m2, LCH_Fe3p_Bottom_Up_m2]...
        = CN_leach(PARAMS, Fe3p_cr, 1, qq, smfull, dz_mm);
    [LCH_TiO4H4_m2, TLCH_TiO4H4_m2, flow_TiO4H4_m2, LCH_TiO4H4_Top_Up_m2, LCH_TiO4H4_Bottom_Down_m2, LCH_TiO4H4_Top_Down_m2, LCH_TiO4H4_Bottom_Up_m2]...
        = CN_leach(PARAMS, TiO4H4_cr, 1, qq, smfull, dz_mm);
    [LCH_Nap_m2, TLCH_Nap_m2, flow_Nap_m2, LCH_Nap_Top_Up_m2, LCH_Nap_Bottom_Down_m2, LCH_Nap_Top_Down_m2, LCH_Nap_Bottom_Up_m2]...
        = CN_leach(PARAMS, Nap_cr, 1, qq, smfull, dz_mm);
    
    % Compute leaching in units of [gr / m^3 / d]
    LCH_Hplus   = LCH_Hplus_m2./(dz_mm/1000);       % [gr/m^3/d]
    dHp_leach1(:,day) = LCH_Hplus;
        LCH_Hplus_Top_Up_m3   = LCH_Hplus_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Hplus_Bottom_Down_m3   = LCH_Hplus_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Hplus_Top_Down_m3   = LCH_Hplus_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Hplus_Bottom_Up_m3   = LCH_Hplus_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
     
    LCH_glucose    = LCH_glucose_m2./(dz_mm/1000);        % [gr/m^3/d]
    dglu_leach1(:,day) = LCH_glucose;
        LCH_glucose_Top_Up_m3   = LCH_glucose_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_glucose_Bottom_Down_m3   = LCH_glucose_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_glucose_Top_Down_m3   = LCH_glucose_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_glucose_Bottom_Up_m3   = LCH_glucose_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        
    LCH_NH3  = LCH_NH3_m2./(dz_mm/1000);      % [gr/m^3/d]
        LCH_NH3_Top_Up_m3   = LCH_NH3_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_NH3_Bottom_Down_m3   = LCH_NH3_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_NH3_Top_Down_m3   = LCH_NH3_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_NH3_Bottom_Up_m3   = LCH_NH3_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        
    LCH_Amm  = LCH_Amm_m2./(dz_mm/1000);      % [gr/m^3/d]
        damm_leach1(:,day) = LCH_Amm;
        LCH_Amm_Top_Up_m3   = LCH_Amm_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Amm_Bottom_Down_m3   = LCH_Amm_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Amm_Top_Down_m3   = LCH_Amm_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Amm_Bottom_Up_m3   = LCH_Amm_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
   
    LCH_Nit     = LCH_Nit_m2./(dz_mm/1000);         % [gr/m^3/d]
        dnit_leach1(:,day) = LCH_Nit;
        LCH_Nit_Top_Up_m3   = LCH_Nit_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Nit_Bottom_Down_m3   = LCH_Nit_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Nit_Top_Down_m3   = LCH_Nit_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Nit_Bottom_Up_m3   = LCH_Nit_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        
    LCH_O2      = LCH_O2_m2./(dz_mm/1000);          % [gr/m^3/d]
        LCH_O2_Top_Up_m3   = LCH_O2_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_O2_Bottom_Down_m3   = LCH_O2_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_O2_Top_Down_m3   = LCH_O2_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_O2_Bottom_Up_m3   = LCH_O2_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        
    LCH_CO2     = LCH_CO2_m2./(dz_mm/1000);         % [gr/m^3/d]
        LCH_CO2_Top_Up_m3   = LCH_CO2_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_CO2_Bottom_Down_m3   = LCH_CO2_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_CO2_Top_Down_m3   = LCH_CO2_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_CO2_Bottom_Up_m3   = LCH_CO2_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        
    LCH_HCO3     = LCH_HCO3_m2./(dz_mm/1000);         % [gr/m^3/d]
        LCH_HCO3_Top_Up_m3   = LCH_HCO3_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_HCO3_Bottom_Down_m3   = LCH_HCO3_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_HCO3_Top_Down_m3   = LCH_HCO3_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_HCO3_Bottom_Up_m3   = LCH_HCO3_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        
    LCH_CO32m     = LCH_CO32m_m2./(dz_mm/1000);         % [gr/m^3/d]
        LCH_CO32m_Top_Up_m3   = LCH_CO32m_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_CO32m_Bottom_Down_m3   = LCH_CO32m_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_CO32m_Top_Down_m3   = LCH_CO32m_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_CO32m_Bottom_Up_m3   = LCH_CO32m_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
                    
    LCH_N2     = LCH_N2_m2./(dz_mm/1000);         % [gr/m^3/d]
        LCH_N2_Top_Up_m3   = LCH_N2_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_N2_Bottom_Down_m3   = LCH_N2_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_N2_Top_Down_m3   = LCH_N2_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_N2_Bottom_Up_m3   = LCH_N2_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        
    LCH_Cbaq     = LCH_Cbaq_m2./(dz_mm/1000);         % [gr/m^3/d]
        LCH_Cbaq_Top_Up_m3   = LCH_Cbaq_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Cbaq_Bottom_Down_m3   = LCH_Cbaq_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Cbaq_Top_Down_m3   = LCH_Cbaq_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Cbaq_Bottom_Up_m3   = LCH_Cbaq_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
    
    LCH_aceticacid     = LCH_aceticacid_m2./(dz_mm/1000);         % [gr/m^3/d]
        LCH_aceticacid_Top_Up_m3   = LCH_aceticacid_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_aceticacid_Bottom_Down_m3   = LCH_aceticacid_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_aceticacid_Top_Down_m3   = LCH_aceticacid_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_aceticacid_Bottom_Up_m3   = LCH_aceticacid_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
    
    LCH_SiO2     = LCH_SiO2_m2./(dz_mm/1000);         % [gr/m^3/d]
        LCH_SiO2_Top_Up_m3   = LCH_SiO2_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_SiO2_Bottom_Down_m3   = LCH_SiO2_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_SiO2_Top_Down_m3   = LCH_SiO2_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_SiO2_Bottom_Up_m3   = LCH_SiO2_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        
    LCH_Al3p     = LCH_Al3p_m2./(dz_mm/1000);         % [gr/m^3/d]
        LCH_Al3p_Top_Up_m3   = LCH_Al3p_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Al3p_Bottom_Down_m3   = LCH_Al3p_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Al3p_Top_Down_m3   = LCH_Al3p_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Al3p_Bottom_Up_m3   = LCH_Al3p_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        
    LCH_Mg2p     = LCH_Mg2p_m2./(dz_mm/1000);         % [gr/m^3/d]
        LCH_Mg2p_Top_Up_m3   = LCH_Mg2p_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Mg2p_Bottom_Down_m3   = LCH_Mg2p_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Mg2p_Top_Down_m3   = LCH_Mg2p_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Mg2p_Bottom_Up_m3   = LCH_Mg2p_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        
    LCH_Ca2p     = LCH_Ca2p_m2./(dz_mm/1000);         % [gr/m^3/d]
        LCH_Ca2p_Top_Up_m3   = LCH_Ca2p_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Ca2p_Bottom_Down_m3   = LCH_Ca2p_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Ca2p_Top_Down_m3   = LCH_Ca2p_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Ca2p_Bottom_Up_m3   = LCH_Ca2p_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        
    LCH_Kp     = LCH_Kp_m2./(dz_mm/1000);         % [gr/m^3/d]
        LCH_Kp_Top_Up_m3   = LCH_Kp_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Kp_Bottom_Down_m3   = LCH_Kp_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Kp_Top_Down_m3   = LCH_Kp_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Kp_Bottom_Up_m3   = LCH_Kp_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        
    LCH_Fe2p     = LCH_Fe2p_m2./(dz_mm/1000);         % [gr/m^3/d]
        LCH_Fe2p_Top_Up_m3   = LCH_Fe2p_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Fe2p_Bottom_Down_m3   = LCH_Fe2p_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Fe2p_Top_Down_m3   = LCH_Fe2p_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Fe2p_Bottom_Up_m3   = LCH_Fe2p_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]        
    
    LCH_Fe3p     = LCH_Fe3p_m2./(dz_mm/1000);         % [gr/m^3/d]
        LCH_Fe3p_Top_Up_m3   = LCH_Fe3p_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Fe3p_Bottom_Down_m3   = LCH_Fe3p_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Fe3p_Top_Down_m3   = LCH_Fe3p_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Fe3p_Bottom_Up_m3   = LCH_Fe3p_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]        

    LCH_TiO4H4     = LCH_TiO4H4_m2./(dz_mm/1000);         % [gr/m^3/d]
        LCH_TiO4H4_Top_Up_m3   = LCH_TiO4H4_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_TiO4H4_Bottom_Down_m3   = LCH_TiO4H4_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_TiO4H4_Top_Down_m3   = LCH_TiO4H4_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_TiO4H4_Bottom_Up_m3   = LCH_TiO4H4_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
                
    LCH_Nap     = LCH_Nap_m2./(dz_mm/1000);         % [gr/m^3/d]
        LCH_Nap_Top_Up_m3   = LCH_Nap_Top_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Nap_Bottom_Down_m3   = LCH_Nap_Bottom_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Nap_Top_Down_m3   = LCH_Nap_Top_Down_m2./(dz_mm/1000);       % [gr/m^3/d]
        LCH_Nap_Bottom_Up_m3   = LCH_Nap_Bottom_Up_m2./(dz_mm/1000);       % [gr/m^3/d]
        
        
%*************************************************************************
%      STATE UPDATES 
%************ *************************************************************
    %%-------------------------------------------------------------------%%
%     CN_solve_negative_N();  

[LCH_Hplus_Top_Up_m3,LCH_Hplus_Bottom_Up_m3,LCH_Hplus,LCH_Hplus_m2] = ...
    fixnegative_solute(Hplus_cr,LCH_Hplus_Top_Up_m3,LCH_Hplus_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_glucose_Top_Up_m3,LCH_glucose_Bottom_Up_m3,LCH_glucose,LCH_glucose_m2] = ...
    fixnegative_solute(glucose_cr,LCH_glucose_Top_Up_m3,LCH_glucose_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_O2_Top_Up_m3,LCH_O2_Bottom_Up_m3,LCH_O2,LCH_O2_m2] = ...
    fixnegative_solute(O2_cr,LCH_O2_Top_Up_m3,LCH_O2_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_CO2_Top_Up_m3,LCH_CO2_Bottom_Up_m3,LCH_CO2,LCH_CO2_m2] = ...
    fixnegative_solute(CO2_gl,LCH_CO2_Top_Up_m3,LCH_CO2_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_HCO3_Top_Up_m3,LCH_HCO3_Bottom_Up_m3,LCH_HCO3,LCH_HCO3_m2] = ...
    fixnegative_solute(HCO3_cr,LCH_HCO3_Top_Up_m3,LCH_HCO3_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_CO32m_Top_Up_m3,LCH_CO32m_Bottom_Up_m3,LCH_CO32m,LCH_CO32m_m2] = ...
    fixnegative_solute(CO32m_cr,LCH_CO32m_Top_Up_m3,LCH_CO32m_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_N2_Top_Up_m3,LCH_N2_Bottom_Up_m3,LCH_N2,LCH_N2_m2] = ...
    fixnegative_solute(N2_cr,LCH_N2_Top_Up_m3,LCH_N2_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_Cbaq_Top_Up_m3,LCH_Cbaq_Bottom_Up_m3,LCH_Cbaq,LCH_Cbaq_m2] = ...
    fixnegative_solute(Cbaq_cr,LCH_Cbaq_Top_Up_m3,LCH_Cbaq_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_aceticacid_Top_Up_m3,LCH_aceticacid_Bottom_Up_m3,LCH_aceticacid,LCH_aceticacid_m2] = ...
    fixnegative_solute(aceticacid_cr,LCH_aceticacid_Top_Up_m3,LCH_aceticacid_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_Amm_Top_Up_m3,LCH_Amm_Bottom_Up_m3,LCH_Amm,LCH_Amm_m2] = ...
    fixnegative_solute(Amm_gl,LCH_Amm_Top_Up_m3,LCH_Amm_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_Nit_Top_Up_m3,LCH_Nit_Bottom_Up_m3,LCH_Nit,LCH_Nit_m2] = ...
    fixnegative_solute(Nit_gl,LCH_Nit_Top_Up_m3,LCH_Nit_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_NH3_Top_Up_m3,LCH_NH3_Bottom_Up_m3,LCH_NH3,LCH_NH3_m2] = ...
    fixnegative_solute(NH3_gl,LCH_NH3_Top_Up_m3,LCH_NH3_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_SiO2_Top_Up_m3,LCH_SiO2_Bottom_Up_m3,LCH_SiO2,LCH_SiO2_m2] = ...
    fixnegative_solute(SiO2_cr,LCH_SiO2_Top_Up_m3,LCH_SiO2_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_Al3p_Top_Up_m3,LCH_Al3p_Bottom_Up_m3,LCH_Al3p,LCH_Al3p_m2] = ...
    fixnegative_solute(Al3p_cr,LCH_Al3p_Top_Up_m3,LCH_Al3p_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_Mg2p_Top_Up_m3,LCH_Mg2p_Bottom_Up_m3,LCH_Mg2p,LCH_Mg2p_m2] = ...
    fixnegative_solute(Mg2p_cr,LCH_Mg2p_Top_Up_m3,LCH_Mg2p_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_Ca2p_Top_Up_m3,LCH_Ca2p_Bottom_Up_m3,LCH_Ca2p,LCH_Ca2p_m2] = ...
    fixnegative_solute(Ca2p_cr,LCH_Ca2p_Top_Up_m3,LCH_Ca2p_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_Kp_Top_Up_m3,LCH_Kp_Bottom_Up_m3,LCH_Kp,LCH_Kp_m2] = ...
    fixnegative_solute(Kp_cr,LCH_Kp_Top_Up_m3,LCH_Kp_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_Fe2p_Top_Up_m3,LCH_Fe2p_Bottom_Up_m3,LCH_Fe2p,LCH_Fe2p_m2] = ...
    fixnegative_solute(Fe2p_cr,LCH_Fe2p_Top_Up_m3,LCH_Fe2p_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_Fe3p_Top_Up_m3,LCH_Fe3p_Bottom_Up_m3,LCH_Fe3p,LCH_Fe3p_m2] = ...
    fixnegative_solute(Fe3p_cr,LCH_Fe3p_Top_Up_m3,LCH_Fe3p_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_TiO4H4_Top_Up_m3,LCH_TiO4H4_Bottom_Up_m3,LCH_TiO4H4,LCH_TiO4H4_m2] = ...
    fixnegative_solute(TiO4H4_cr,LCH_TiO4H4_Top_Up_m3,LCH_TiO4H4_Bottom_Down_m3,dz_mm, nl_soil);
[LCH_Nap_Top_Up_m3,LCH_Nap_Bottom_Up_m3,LCH_Nap,LCH_Nap_m2] = ...
    fixnegative_solute(Nap_cr,LCH_Nap_Top_Up_m3,LCH_Nap_Bottom_Down_m3,dz_mm, nl_soil);

% CN_solve_negative_N_SRM(); 
      damm_leach2(:,day) = LCH_Amm;
      dnit_leach2(:,day) = LCH_Nit;
      dglu_leach2(:,day) = LCH_glucose;
      dHp_leach2(:,day) = LCH_Hplus;
  
    %%-------------------------------------------------------------------%% 
    % Update all concentrations
    omega = 10^-14;
%     H+
        dHplus_dt = rainadd_Hp(:,day) - LCH_Hplus;%;
        Hplus_new = Hplus_cr + dHplus_dt;

%     C6H12O6
        dglucose_dt =  - LCH_glucose;%;
        glucose_new = glucose_cr + dglucose_dt; %what is left in soil; cexg is amount exuded
        glucose_new(abs(glucose_new)<=omega)=0;  

%     NH3(aq) 
        dNH3_dt = - LCH_NH3;%;
        NH3_new = NH3_gl + dNH3_dt;
        NH3_new(abs(NH3_new)<=omega)=0; 
               
%     NH4+
        dAmm_dt =  -LCH_Amm;%; 
        Amm_new = Amm_gl + dAmm_dt;
        % To avoid non-physical negative values
        Amm_new(Amm_new<0) = 1e-10;

%     NO3-
        dNit_dt = rainadd_NO3m(:,day) -LCH_Nit;% ;
        Nit_new = Nit_gl + dNit_dt;
        % To avoid non-physical negative values
        Nit_new(Nit_new<0) = 1e-10;
        Nit_new(abs(Nit_new)<=omega)=0;  

%     O2(aq)
        dO2_dt = rainadd_O2(:,day) - LCH_O2;%;
        O2_new = O2_cr + dO2_dt;
        O2_new(abs(O2_new)<=omega)=1e-11;    

%     CO2(aq)
        dCO2_dt = -LCH_CO2;%;
        CO2_new = CO2_cr + dCO2_dt;
        CO2_new(abs(CO2_new)<=omega)=0;  

%      HCO3
        dHCO3_dt = rainadd_HCO3m(:,day) - LCH_HCO3;%;
        HCO3_new = HCO3_cr + dHCO3_dt;
        HCO3_new(abs(HCO3_new)<=omega)=0;  

%      CO32m
        dCO32m_dt = -LCH_CO32m;%;
        CO32m_new = CO32m_cr + dCO32m_dt;
        CO32m_new(abs(CO32m_new)<=omega)=0; 
        
 %     N2
        dN2_dt = -LCH_N2;%;
        N2_new = N2_cr + dN2_dt;
        N2_new(abs(N2_new)<=omega)=0;  
        
%     Cb_aq
        dCbaq_dt = -LCH_Cbaq;%;
        Cbaq_new = Cbaq_cr + dCbaq_dt;
        Cbaq_new(abs(Cbaq_new)<=omega)=0;  

%     Acetic Acid
        daceticacid_dt = -LCH_aceticacid;%;
        aceticacid_new = aceticacid_cr + daceticacid_dt;
        toolow = find(aceticacid_new<0.5);
        aceticacid_new(aceticacid_new<0.5)=1;  
        aceticacidoverride(toolow) = 1;        
                
%     SiO2
        dSiO2_dt = -LCH_SiO2;%;
        SiO2_new = SiO2_cr + dSiO2_dt;
        SiO2_new(abs(SiO2_new)<=omega)=0;          
                
%     Al3p
        dAl3p_dt = -LCH_Al3p;%;
        Al3p_new = Al3p_cr + dAl3p_dt;
        Al3p_new(abs(Al3p_new)<=omega)=0;          
                
%     Mg2p
        dMg2p_dt = rainadd_Mg2p(:,day) - LCH_Mg2p;%;
        Mg2p_new = Mg2p_cr + dMg2p_dt;
        Mg2p_new(abs(Mg2p_new)<=omega)=0;          
                
%     Ca2p
        dCa2p_dt = rainadd_Ca2p(:,day)-LCH_Ca2p;%;
        Ca2p_new = Ca2p_cr + dCa2p_dt;
        Ca2p_new(abs(Ca2p_new)<=omega)=0;          
                
%     Kp
        dKp_dt = rainadd_Kp(:,day)-LCH_Kp;%;
        Kp_new = Kp_cr + dKp_dt;
        Kp_new(abs(Kp_new)<=omega)=0;          
                
%     Fe2p
        dFe2p_dt = -LCH_Fe2p;%;
        Fe2p_new = Fe2p_cr + dFe2p_dt;
        Fe2p_new(abs(Fe2p_new)<=omega)=0;          

%     Fe3p
        dFe3p_dt = -LCH_Fe3p;%;
        Fe3p_new = Fe3p_cr + dFe3p_dt;
        Fe3p_new(abs(Fe3p_new)<=omega)=0; 
        
%     TiO4H4
        dTiO4H4_dt = -LCH_TiO4H4;%;
        TiO4H4_new = TiO4H4_cr + dTiO4H4_dt;
        TiO4H4_new(abs(TiO4H4_new)<=omega)=0;          
                
%     Nap
        dNap_dt = rainadd_Nap(:,day)-LCH_Nap;%;
        Nap_new = Nap_cr + dNap_dt;
        Nap_new(abs(Nap_new)<=omega)=0;         
        

    Cl_new = Cl + dCl_dt*dtime/86400;
    
    if SWITCHES.CN_type
        Ch_new = Ch + dCh_dt*dtime/86400;
    else
        Ch_new = NaN;
    end
    
  
    Cb_new = Cb_cr + dCb_dt;          %[gC/m3]
   
    indcb = Cb_new<0; 
    Cb_new(indcb) = 0;    
    
% %*************************************************************************
% %      Mass balance of C and N
% %*************************************************************************

    if SWITCHES.CN.Bioturbation
        % CHECK MASS BALANCE OF C
        dz = [VERTSTRUC.dzsmm/1000]; %[m] grid depth
        Cinput = sum(ADD.*dz);       %[grC/m2/d] Input from Litter and Bioturbation
        CBoutput = sum(OUT.*dz);     %[grC/m2/d] Output from Bioturbation
        if SWITCHES.CN.Bioturbation_new == 1
            if SWITCHES.CN_type == 1
                Cinput = sum(ADD.*dz + VARIABLES.Cb_bio_in.*dz + VARIABLES.Ch_bio_in.*dz);
                CBoutput = sum(OUT.*dz + VARIABLES.Cb_bio_out.*dz + VARIABLES.Ch_bio_out.*dz);
            elseif SWITCHES.CN_type == 0
                Cinput = sum(ADD.*dz + VARIABLES.Cb_bio_in.*dz);
                CBoutput = sum(OUT.*dz + VARIABLES.Cb_bio_out.*dz);
            end
        end
        DECout = sum(DECl.*dz*rr);   %[grC/m2/d] Respiration
        dCb_m2 = sum(dCb_dt.*dz);    %[grC/m2/d] change in Microbial Biomass 
        dCl_m2 = sum(dCl_dt.*dz);    %[grC/m2/d] change in organic matter
        mberrorC = Cinput - CBoutput - DECout - dCb_m2 - dCl_m2; %[grC/m2/d]

%         DECout = sum(DECl.*dz*rr);   %[grC/m2/d] Respiration
                DECout = sum(DECl.*dz);   %[grC/m2/d] Respiration
        dCb_m2 = sum(dCb_dt.*dz);    %[grC/m2/d] change in Microbial Biomass 
        dCl_m2 = sum(dCl_dt.*dz);    %[grC/m2/d] change in organic matter
        mberrorC = Cinput - CBoutput - DECout - dCb_m2 - dCl_m2; %[grC/m2/d]
        
        % CHECK MASS BALANCE OF N
        Ninput = sum(ADD./CNa.*dz);      %[gr/m2/d]
        NBoutput = sum(OUT./CNo.*dz);    %[gr/m2/d]  
        if SWITCHES.CN.Bioturbation_new == 1
            if SWITCHES.CN_type == 1
                Ninput = sum(ADD./CNa.*dz + VARIABLES.Cb_bio_in./CNb.*dz + VARIABLES.Ch_bio_in./VARIABLES.Ch_CN_bio_in.*dz);
                NBoutput = sum(OUT./CNo.*dz + VARIABLES.Cb_bio_out./CNb.*dz + VARIABLES.Ch_bio_out./VARIABLES.Ch_CN_bio_out.*dz);
            elseif SWITCHES.CN_type == 0
                Ninput = sum(ADD./CNa.*dz + VARIABLES.Cb_bio_in./CNb.*dz);
                NBoutput = sum(OUT./CNo.*dz + VARIABLES.Cb_bio_out./CNb.*dz);                
            end
        end
        dNb_dt = dCb_dt./113.1146.*14.01;
        dAmm_m2 = sum(dAmm_dt.*dz);      %[gr/m2/d]
        dNit_m2 = sum(dNit_dt.*dz);      %[gr/m2/d]
        dNH3_m2 = sum(dNH3_dt.*dz);      %[gr/m2/d]
        dNb_m2 = sum(dNb_dt.*dz);        %[gr/m2/d]
                dNl_m2 = sum(dz*0);        %[gr/m2/d]
        inputN = Ninput - NBoutput  ...
            +sum(Depo_Nit_m2+Depo_Amm_m2)+ sum(Urea); %[gr/m2/d] 
        outputN = dNb_m2 + dNl_m2 + dAmm_m2 + dNit_m2 + dNH3_m2; %[gr/m2/d]
        mberrorN = inputN - outputN; %[gr/m2/d]
        
    
    % Bioturbation should happen all the time!!
    else
        % CHECK MASS BALANCE OF C
        Linput = sum(ADD.*(dz_mm/1000)); %[grC/m2/d]
        DECout = sum(DECl.*(dz_mm/1000)*rr);       %[grC/m2/d]
        dCb_m2 = sum(dCb_dt.*(dz_mm/1000));   %[grC/m2/d]
        dCl_m2 = sum(dCl_dt.*(dz_mm/1000));   %[grC/m2/d]
        mberrorC = Linput - DECout - dCb_m2 - dCl_m2; %[grC/m2/d]
        
        % CHECK MASS BALANCE OF N
        Linput = sum(ADD./CNa.*(dz_mm/1000)); %[gr/m2/d]
        dAmm_m2 = sum(dAmm_dt.*(dz_mm/1000)); %[gr/m2/d]
        dNit_m2 = sum(dNit_dt.*(dz_mm/1000)); %[gr/m2/d]
        dNb_m2 = sum(dNb_dt.*(dz_mm/1000));   %[gr/m2/d]
        dNl_m2 = sum(dNl_dt.*(dz_mm/1000));   %[gr/m2/d]
        inputN = Linput - LCH_N_m2;   %- UP_N_m2     %[gr/m2/d]
        outputN = dAmm_m2 + dNit_m2 + dNb_m2 + dNl_m2; %[gr/m2/d]
        mberrorN = inputN - outputN; %[gr/m2/d]
    end

       
%*************************************************************************
%       STORE IN STRUCTURE 
%*************************************************************************
    VARIABLES.Cl=Cl_new; 
    VARIABLES.Ch=Ch_new;    
    VARIABLES.Cb=Cb_new;
    
    VARIABLES.Nit=Nit_new;
    VARIABLES.NH3=NH3_new;
    VARIABLES.Amm=Amm_new;
    
    VARIABLES.pH = pH_cr;
      
    VARIABLES.SUP_amm=SUP_amm;
    VARIABLES.SUP_amm_all=SUP_amm_all;
    VARIABLES.SUP_amm_all_m2=SUP_amm_all_m2;    
    VARIABLES.SUP_nit=SUP_nit;
    VARIABLES.SUP_nit_all=SUP_nit_all;
    VARIABLES.SUP_nit_all_m2=SUP_nit_all_m2;    
    VARIABLES.Tot_UP_N_m2=Tot_UP_N_m2;
    VARIABLES.Tot_UP_Amm=Tot_UP_Amm;
    VARIABLES.Tot_UP_Nit=Tot_UP_Nit;
    
    VARIABLES.LCH_Amm=LCH_Amm; 
    VARIABLES.LCH_Nit=LCH_Nit;
    VARIABLES.TLCH_Amm_m2=TLCH_Amm_m2; 
    VARIABLES.TLCH_Nit_m2=TLCH_Nit_m2;

    VARIABLES.PHI=PHI;
    VARIABLES.phi=phi; 
    VARIABLES.fSd=fSd;
    VARIABLES.fTd=fTd;
    
    % Crunch Variables
    VARIABLES.CO2 = CO2_new; % mol/kg 
    VARIABLES.HCO3 = HCO3_new; % mol/kg 
    VARIABLES.CO32m = CO32m_new; % mol/kg 
    VARIABLES.N2 = N2_new; % mol/kg 
    VARIABLES.O2 = O2_new; % mol/kg 
    VARIABLES.Hplus = Hplus_new; % mol/kg
    VARIABLES.glucose = glucose_new; 
    VARIABLES.Cbaq = Cbaq_new; 
    VARIABLES.SiO2 = SiO2_new; 
    VARIABLES.Al3p = Al3p_new; 
    VARIABLES.Mg2p = Mg2p_new; 
    VARIABLES.Ca2p = Ca2p_new; 
    VARIABLES.Kp = Kp_new; 
    VARIABLES.Fe2p = Fe2p_new; 
    VARIABLES.Fe3p = Fe3p_new; 
    VARIABLES.TiO4H4 = TiO4H4_new; 
    VARIABLES.Nap = Nap_new;    
    VARIABLES.Albite = Albite_cr;
    VARIABLES.Anatase = Anatase_cr;
    VARIABLES.Calcite = Calcite_cr;
    VARIABLES.Chamosite = Chamosite_cr;
    VARIABLES.Anthophyllite = Anthophyllite_cr;
    VARIABLES.Kaolinite = Kaolinite_cr;
    VARIABLES.Microcline = Microcline_cr;  
    VARIABLES.Muscovite = Muscovite_cr;
    VARIABLES.Quartz = Quartz_cr;
    VARIABLES.Goethite = Goethite_cr;
    VARIABLES.Lepidocrocite= Lepidocrocite_cr;
    VARIABLES.Epidote = Epidote_cr;
    VARIABLES.Magnetite = Magnetite_cr;

    
%*************************************************************************
%       Save the necessry variables 
%*************************************************************************
    % Litter Input
    ADD_net_print(:,day)=ADD_net;                % [g/m3/day]
    
    % Organic Carbon
    Cl_print(:,day)=Cl_new;                  % [g/m3/day]
    Cb_print(:,day)=Cb_new;                  % [g/m3/day]
    Ch_print(:,day)=Ch_new;                  % [g/m3/day]

    
    % Organic Nitrogen
    Nb_print(:,day)=Cb_new./CN_biomass;      % [g/m3/day]

    % Inorganic Nitrogen  
    Nit_print(:,day)=Nit_new;                % [g/m3/day]
    NH3_print(:,day)=NH3_new;                % [g/m3/day]
    Amm_print(:,day)=Amm_new;                % [g/m3/day]

    
    % Flux save
    LCH_Hplus_Bottom_Down_m2_print(:,day) = LCH_Hplus_Bottom_Up_m2;
    LCH_Hplus_Bottom_Up_m2_print(:,day) = LCH_Hplus_Bottom_Up_m2;
	LCH_Hplus(:,day) = LCH_Hplus;
    TLCH_Hplus_m2_print(:,day)=TLCH_Hplus_m2;
    
	LCH_glucose_Bottom_Down_m2_print(:,day) = LCH_glucose_Bottom_Up_m2;
    LCH_glucose_Bottom_Up_m2_print(:,day) = LCH_glucose_Bottom_Up_m2;
	TLCH_glucose_m2_print(:,day)=TLCH_glucose_m2;

	LCH_NH3_Bottom_Down_m2_print(:,day) = LCH_NH3_Bottom_Up_m2;
    LCH_NH3_Bottom_Up_m2_print(:,day) = LCH_NH3_Bottom_Up_m2;
        TLCH_NH3_m2_print(:,day)=TLCH_NH3_m2;

    LCH_Amm_Bottom_Down_m2_print(:,day) = LCH_Amm_Bottom_Up_m2;
    LCH_Amm_Bottom_Up_m2_print(:,day) = LCH_Amm_Bottom_Up_m2;
	    TLCH_Amm_m2_print(:,day)=TLCH_Amm_m2;

	LCH_Nit_Bottom_Down_m2_print(:,day) = LCH_Nit_Bottom_Up_m2;
    LCH_Nit_Bottom_Up_m2_print(:,day) = LCH_Nit_Bottom_Up_m2;
	    TLCH_Nit_m2_print(:,day)=TLCH_Nit_m2;

	LCH_O2_Bottom_Down_m2_print(:,day) = LCH_O2_Bottom_Up_m2;
    LCH_O2_Bottom_Up_m2_print(:,day) = LCH_O2_Bottom_Up_m2;
	    TLCH_O2_m2_print(:,day)=TLCH_O2_m2;

	LCH_CO2_Bottom_Down_m2_print(:,day) = LCH_CO2_Bottom_Up_m2;
    LCH_CO2_Bottom_Up_m2_print(:,day) = LCH_CO2_Bottom_Up_m2;
        TLCH_CO2_m2_print(:,day)=TLCH_CO2_m2;

    LCH_CO32m_Bottom_Down_m2_print(:,day) = LCH_CO32m_Bottom_Up_m2;
    LCH_CO32m_Bottom_Up_m2_print(:,day) = LCH_CO32m_Bottom_Up_m2;
        TLCH_CO32m_m2_print(:,day)=TLCH_CO32m_m2;

    LCH_HCO3_Bottom_Down_m2_print(:,day) = LCH_HCO3_Bottom_Up_m2;
    LCH_HCO3_Bottom_Up_m2_print(:,day) = LCH_HCO3_Bottom_Up_m2;
        TLCH_HCO3_m2_print(:,day)=TLCH_HCO3_m2;

    LCH_N2_Bottom_Down_m2_print(:,day) = LCH_N2_Bottom_Up_m2;
    LCH_N2_Bottom_Up_m2_print(:,day) = LCH_N2_Bottom_Up_m2;
    TLCH_N2_m2_print(:,day)=TLCH_N2_m2;

    LCH_Cbaq_Bottom_Down_m2_print(:,day) = LCH_Cbaq_Bottom_Up_m2;
    LCH_Cbaq_Bottom_Up_m2_print(:,day) = LCH_Cbaq_Bottom_Up_m2;
        TLCH_Cbaq_m2_print(:,day)=TLCH_Cbaq_m2;

    LCH_Acetic_acid_Bottom_Down_m2_print(:,day) = LCH_aceticacid_Bottom_Up_m2;
    LCH_Acetic_acid_Bottom_Up_m2_print(:,day) = LCH_aceticacid_Bottom_Up_m2;
        TLCH_aceticacid_m2_print(:,day)=TLCH_aceticacid_m2;

    LCH_NH3_Bottom_Down_m2_print(:,day) = LCH_NH3_Bottom_Down_m2;
    LCH_NH3_Bottom_Up_m2_print(:,day) = LCH_NH3_Bottom_Up_m2;
        TLCH_NH3_m2_print(:,day)=TLCH_NH3_m2;

    TLCH_SiO2_m2_print(:,day)=TLCH_SiO2_m2;
    TLCH_Al3p_m2_print(:,day)=TLCH_Al3p_m2;
    TLCH_Mg2p_m2_print(:,day)=TLCH_Mg2p_m2;
    TLCH_Ca2p_m2_print(:,day)=TLCH_Ca2p_m2;
    TLCH_Kp_m2_print(:,day)=TLCH_Kp_m2;
    TLCH_Fe2p_m2_print(:,day)=TLCH_Fe2p_m2;
    TLCH_Fe3p_m2_print(:,day)=TLCH_Fe3p_m2;
    TLCH_TiO4H4_m2_print(:,day)=TLCH_TiO4H4_m2;
    TLCH_Nap_m2_print(:,day)=TLCH_Nap_m2;
    
    % Nitrate Flux
    LCH_Nit_print(:,day)=LCH_Nit;            % [g/m3/day]
     SUP_nit_print(:,day)=SUP_nit;
    Depo_Nit_m3_print(:,day)=Depo_Nit_m3;
    
    % Ammonium Flux
    LCH_amm_print(:,day)=LCH_Amm;            % [g/m3/day] 
   LCH_amm_m2_print(:,day)=LCH_Amm_m2;
   TLCH_amm_m2_print(:,day)=TLCH_Amm_m2;
    NH4_Urea_m3_print(:,day)=NH4_Urea_m3;
    Ins_Vol_m3_print(:,day)=Ins_Vol_m3;
     SUP_amm_print(:,day)=SUP_amm;
    Depo_Amm_m3_print(:,day)=Depo_Amm_m3;
    Fert_Amm_m3_print(:,day)=Fert_Amm_m3;
    
    
    % Ammonia Flux
    NH3_Urea_m3_print(:,day)=NH3_Urea_m3;
    LCH_NH3_print(:,day)=LCH_NH3;            % [g/m3/day]
   LCH_NH3_m2_print(:,day)=LCH_NH3_m2;
   TLCH_NH3_m2_print(:,day)=TLCH_NH3_m2;
    Volat_print(:,day)=Volat;
    Ins_Vol_print(:,day)=Ins_Vol_m3;
    
   
CO2gaslost_print(:,day) = respCgas + resp2gas;  % total respired C from crunch/exudates and decomposition

% Crunch Extras
    CO2gas_print(:,day) = CO2gas_cr;  % this will change when diffusion incorporated
    O2gas_print(:,day) = O2gas_cr;  % this will change when diffusion incorporated
    CO2_print(:,day) = CO2_new; % mol/kg 
    N2_print(:,day) = N2_new; % mol/kg 
    O2_print(:,day) = O2_new; % mol/kg 
    Hplus_print(:,day)=Hplus_new; % mol/kg
    pH_print(:,day)=pH_cr;
    glucose_print(:,day) = glucose_new; % mol/kg
    Cbaq_print(:,day) = Cbaq_new; % mol/kg
    HCO3_print(:,day) = HCO3_new; % mol/kg
    CO32m_print(:,day) = CO32m_new; % mol/kg
    aceticacid_print(:,day) = aceticacid_new; % mol/kg
    aceticacidoverride_print(:,day) = aceticacidoverride;
    SiO2_print(:,day) = SiO2_new; % mol/kg
    Al3p_print(:,day) = Al3p_new; % mol/kg
    Mg2p_print(:,day) = Mg2p_new; % mol/kg
    Ca2p_print(:,day) = Ca2p_new; % mol/kg
    Kp_print(:,day) = Kp_new; % mol/kg
    Fe2p_print(:,day) = Fe2p_new; % mol/kg
    Fe3p_print(:,day) = Fe3p_new; % mol/kg
    TiO4H4_print(:,day) = TiO4H4_new; % mol/kg
    Nap_print(:,day) = Nap_new; % mol/kg
    Albite_print(:,day) = Albite_cr;
    Anatase_print(:,day) = Anatase_cr;
    Calcite_print(:,day) = Calcite_cr;
    Chamosite_print(:,day) = Chamosite_cr;
    Anthophyllite_print(:,day) = Anthophyllite_cr;
    Kaolinite_print(:,day) = Kaolinite_cr;
    Microcline_print(:,day) = Microcline_cr;  
    Muscovite_print(:,day) = Muscovite_cr;
    Quartz_print(:,day) = Quartz_cr;
    Goethite_print(:,day) = Goethite_cr;
    Lepidocrocite_print(:,day) = Lepidocrocite_cr;
    Epidote_print(:,day) = Epidote_cr;
    Magnetite_print(:,day) = Magnetite_cr;
    
    SUP_nit_print(:,day)=SUP_nit;            % [g/m3/day]
    SUP_amm_print(:,day)=SUP_amm;            % [g/m3/day]
%     
    SUP_amm_active_print(:,day)=SUP_amm_active; % [g/m3/day]
    SUP_nit_active_print(:,day)=SUP_nit_active; % [g/m3/day]
    SUP_amm_passive_print(:,day)=SUP_amm_passive; % [g/m3/day]
    SUP_nit_passive_print(:,day)=SUP_nit_passive; % [g/m3/day]
%     
    
   Tot_UP_Amm_print(:,day)=Tot_UP_Amm;      % [g/m3/day]
   Tot_UP_Nit_print(:,day)=Tot_UP_Nit;      % [g/m3/day]
    
    if SWITCHES.CN.exudation == 1
        %Exudation
        Gtot_print(:,day)=Gtot; % [g/m3/day]
        cexg_print(:,day)=cexg; % 
        cexf_print(:,day)=cexf; % 
        fflav_print(:,day)=fflav; % 
        Cbglu_print(:,day) = Cbglu;
        Bp_print(:,day) = Bp;
        Br_print(:,day) = Br;

    end
%%-----------------------------------------------------------------------%%
    if save_whole_cycle == 1

        PHI_print(:,day)=PHI;
        phi_print(:,day)=phi;
        fSd_print(:,day)=fSd;
        fTd_print(:,day)=fTd;

	Cbiorate_print(day)=VARIABLES.Cbiorate;
        Nbiorate_print(day)=VARIABLES.Nbiorate;
        bioNerror_print(day)=VARIABLES.bioNerror;
        bioCerror_print(day)=VARIABLES.bioCerror;

        Cl_bio_change_print(:,day)=VARIABLES.Cl_bio_change;
        Cin_m3_print(:,day)=VARIABLES.Cl_bio_in;
        Cout_m3_print(:,day)=VARIABLES.Cl_bio_out;

        flow_Nit_m2_print(:,day)=flow_Nit_m2;
        flow_Amm_m2_print(:,day)=flow_Amm_m2;


    end
end

%%-----------------------------------------------------------------------%%

disp('Simulation End')
toc


