function [Fert_Amm_m3, Fert_Nit_m3, Ins_Vol_m2, NH3_Urea_m2, NH4_Urea_m2, VARIABLES, SWITCHES]...
    = CN_fertilizer(VARIABLES, VERTSTRUC, SWITCHES, day, IN_Nit_Fertilizer,...
    IN_Amm_Fertilizer, IN_Organic_Fertilizer, nl_soil, PPT_mm, Wind_ms, BulkDensity,...
    sm, Ts, smp, WFPS, Clay, pH) 

%==========================================================================  
% Fertilizer application
% - Urea hydrolysis: 
% Grant et al., 2010: Changes in net ecosystem productivity and greenhouse 
% gas exchange with fertilization of Douglas fir: Mathematical modeling in ecosys
% Vlek and Carter, 1983: The effect of soil environment and fertilizer
% modifications on the rate of urea hydrolysis
% Lai et al., 1993: Kinetics of urea hydrolysis in wheat residue
% Metivier et al., 2009: Using the ecosys mathematical model to simulate
% temporal variability of nitrous oxide emissions from a fertilized
% agricultural soil.
% Grant, 2004: Modeling topographic effects on net ecosystem productivity
% of boreal black spruce forests.
% Grant and Rochette, 1994: Soil microbial respiration at different water
% potentials and temperatures: Theory and mathematical modeling
% Moyo et al., 1989: Temperature effects on soil urease activity
% Grant et al., 1993: Simulation of carbon and nitrogen transformations in
% soil: mineralization
%
% - Volatilization
% Sherlock and Goh, 1984: Dynamics of ammonia volatilization from simulated
% urine patches and aqueous urea applied to pasture. II. Theoretical
% derivation of a simplified model.
% Waldrip et al., 2013: Estimation of ammonia emissons from beef cattle
% feedyards using the process-based model manure-DNDC
% Bouwmeester and Vlek, 1980: Rate control of ammonia volatilization from
% rice paddies.
% Li, 2000: Modeling trace gas emissions from agricultural ecosystems
% Li et al., 1992: A model of nitrous oxide evolution from soil driven by
% rainfall events: 1. Model structure and sensitivity
% Sommer et al., 2004: Ammonia emisson from mineral fertilizers and
% fertilized crops
%
% - Modeling approach
% Overdahl and Rehm, 1990: Using anhydrous ammonia in minnesota
% Li et al., 1994: Modeling carbon biogeochemistry in agricultural soils
% Preez and Burger, 1988: Movement of ammonia plus ammonium from nitrogen
% fertilizers band placed in alkaline soils
%
% Written by Dongkook Woo, UIUC, 2014
% All rights reserved!

%------------------------- Input Variables --------------------------------
% day                   % [d] Time step, here day.
% IN_Nit_Fertilizer     % [gN/m2] Nitrate fertilizer
% IN_Amm_Fertilizer     % [gN/m2] Ammonium fertilizer
% IN_Organic_Fertilizer % [gN/m2] Urea fertilizer
% Cb                    % [gC/m3] Carbon concentration in biomass pool 
% dz_mm                 % [mm] Grid depth 
% nl_soil               % number of soil layer
% Urea                  % [gN/m2] Urea
% PPT_mm                % [mm] Precipitation
% Wind_ms               % [m/s] Wind speed
% Hydrolysis_on         % Hydrolysis switch
% BulkDensity           % [g/m3] Bulk density
% sm                    % [-] soil moisture
% Ts                    % ['C] Soil temperature
% smp                   % [MPa] Soil matric potential
% NH3                   % [gN/m3] Ammonia
% WFPS                  % [-] Water Filled Porosity 
% Clay                  % [-] Clay amount
Cb = VARIABLES.Cb;   
dz_mm=VERTSTRUC.dzsmm;
Hydrolysis_on=SWITCHES.CN.Hydrolysis_on; 
Urea=VARIABLES.Urea;

%------------------------- Parameters -------------------------------------
Dsi=1;                   % [gC/g micr. C/h] Rate constant for hydrolysis for protein
% http://www.ingredients101.com/urea.htm
% Urea is a synthetic compound with a high level of nitrogen. It is not a 
% natural protein. In high-energy rations, it can serve as a source of 
% nitrogen for the synthesis of protein by rumen microbes. This is why it 
% is said to have an equivalent crude protein value from nonprotein nitrogen.
KmD=75;                  % [gC/Mg] Michaelis-Menten constant for hydrolysis
KiD=25;                  % [gC/m3-liquid] Inhibition constant for hydrolysis, typo (*(1/1000000)*BulkDensity(1)/sm(1);)

A=17.1124;               % [-] Parameter selected such that ftg=1 at TsK=30'C  
Ha=57500;                % [J/mol] Energy of activation
R=8.3295;                % [J/mol/K] Universal gas constant
Hdl=192500;              % [J/mol] Energy of low temperature deactivation
E=710;                   % [J/mol/K] Change in entropy
Hdh=222500;              % [J/mol] Energy of high temperature deactivation

Time_adj=24/1;           % [h] Time step[day]/hour

%------------------------- Output Variables -------------------------------
% Fert_Amm_m3            % [gN/m3/d] Ammonium fertilizer input   
% Fert_Nit_m3            % [gN/m3/d] Nitrate fertilizer input  
% VARIABLES.Urea         % [gN/m2] Remaining urea
% SWITCHES
% VARIABLES.NH3          % [gN/m2/d] NH3-N 
% Ins_Vol                % [gN/m2/d] Instaneous volatilization

%==========================================================================    
% Initialization
Fert_Nit_m3=zeros(nl_soil,1);
Fert_Amm_m3=zeros(nl_soil,1);

%==========================================================================    
% Ammonium fertilizer [gN/m2] -> [gN/m3]
Fert_Amm_m3(1)=IN_Amm_Fertilizer(day)/(1)/(dz_mm(1)./1000);
%Fert_Amm_m3(2)=IN_Amm_Fertilizer(day)/(3)/(dz_mm(2)./1000);
%Fert_Amm_m3(3)=IN_Amm_Fertilizer(day)/(3)/(dz_mm(3)./1000);

%==========================================================================    
% Nitrate fertilizer [gN/m2] -> [gN/m3]
Fert_Nit_m3(1)=IN_Nit_Fertilizer(day)/(1)/(dz_mm(1)./1000);
%Fert_Nit_m3(2)=IN_Nit_Fertilizer(day)/(3)/(dz_mm(2)./1000);
%Fert_Nit_m3(3)=IN_Nit_Fertilizer(day)/(3)/(dz_mm(3)./1000);

%==========================================================================    
% Urea hydrolysis 
Urea=Urea+IN_Organic_Fertilizer(day); % [gN/m2]
Urea_m3=Urea./(dz_mm(1)./1000); % [gN/m3]

% Check precipitation occurance
if Urea>0
    if PPT_mm(day) > 0
        Hydrolysis_on=1; % on
    end
elseif Urea <= 0
    Hydrolysis_on=0; % off
end

% Estimate urea hydrolysis 
% CO(NH2)2,l+H2O->CO2,l+2NH3,l

% http://www.convertunits.com/molarmass/CO(NH2)2
% Molecular weight calculation: 12.0107 + 15.9994 + (14.0067 + 1.00794*2)*2

%Element  Symbol Atomic  Mass # of Atoms Mass Percent
%Hydrogen H      1.00794      4          6.713%
%Carbon   C      12.0107      1	         19.999%
%Nitrogen N      14.0067      2          46.646%
%Oxygen   O      15.9994      1          26.641% 

% Molar mass of CO(NH2)2 = 60.05526 g/mol
Ins_Vol_m2=0;
NH3_Urea_m2=zeros(nl_soil,1);
NH4_Urea_m2=zeros(nl_soil,1);
D_N_d=0;
if Hydrolysis_on == 1;
    Urea_m3_C=Urea_m3*(19.999/46.646); %[gC/m3]
    Urea_g_C=Urea_m3_C./BulkDensity(1); % [gC/g-soil]
    Urea_Mg_C=Urea_g_C*(1000000/1); % [gC/Mg-soil]
    S=Urea_Mg_C;
    
    Cb_liquid=Cb(1)./sm(1); % [gC/m3-liquid]
    M_l=Cb_liquid;
    M=Cb(1).*(dz_mm(1)./1000); % [gC/m2]
    
    DD=(Dsi*S)./(S+KmD*(1+(M_l/KiD))); % [gC/g micro.C/h]
    
    % Temperature dependent function
    TsK=Ts(1)+273.15; % ['C] -> [K]
    
    ftg=TsK*((exp(A-Ha/(R*TsK))) / (1+exp((Hdl-E*TsK)/(R*TsK))+exp((E*TsK-Hdh)/(R*TsK))));
    indGT1=ftg>1;
    ftg(indGT1)=1;
    indLT0=ftg<0;
    ftg(indLT0)=0;
    
    % Moisture dependent function
    skp=smp(1)*(1000/1); 
    if skp > 20
        fsg=-0.11*log(skp)+1.3269;
    else
        fsg=0.0173*log(skp)+0.9475;        
    end
    indGT1=fsg>1;
    fsg(indGT1)=1;
    indLT0=fsg<0;
    fsg(indLT0)=0;

    % The rate at which C is hydrolyzed [gC/m2/h]
    D=DD*M*ftg*fsg;
    
    % Scaling up.
    D_d=D*Time_adj; %[gC/m2/d]
    
    % Carbon to Nitrogen
    % CO2,l+2NH3,l
    %Element  Symbol Atomic  Mass # of Atoms Mass Percent
    %Hydrogen H      1.00794      6          7.746%
    %Carbon   C      12.0107      1	         15.384%
    %Nitrogen N      14.0067      2          35.882%
    %Oxygen   O      15.9994      2          40.987% 
    D_N_d=D_d*(35.882/15.384); % NH3-N [gN/m2/d]
    if D_N_d > Urea
        D_N_d=Urea;
    end
    
    %---------------------------------------------------------------------
    % Instaneous volatilization
    [Ins_Vol_m2] = CN_volatilization(D_N_d, Ts, WFPS, Clay);
    D_N_d = D_N_d - Ins_Vol_m2;
    %---------------------------------------------------------------------
    % Check
    %Urea+NH3_Update+Ins_Vol % This value should be the same as pre-urea
end

%--------------------------------------------------------------------------
% Instaneous exchange between NH3 and NH4
CN_NH4m_fraction ();

% As can be shown, there are only one N in NH4 and NH3. Thus, we can use
% the fraction of NH4m directly
NH3_NH4=D_N_d;

NH4_Urea=NH3_NH4*fraction_NH4m(1); % [gN/m2]
NH3_Urea=NH3_NH4*(1-fraction_NH4m(1)); % [gN/m2]

% Check the mass balance
%checkMbarrance=NH4_Urea+NH3_Urea-NH3_NH4;


%--------------------------------------------------------------------------
% Update
VARIABLES.Urea=Urea-D_N_d;
SWITCHES.CN.Hydrolysis_on=Hydrolysis_on;

NH3_Urea_m2(1)=NH3_Urea; % NH3-N [g/m2/d]
%NH3_Urea_m2(2)=NH3_Urea/3; % NH3-N [g/m2/d]
%NH3_Urea_m2(3)=NH3_Urea/3; % NH3-N [g/m2/d]
NH4_Urea_m2(1)=NH4_Urea;
%NH4_Urea_m2(2)=NH4_Urea/3;
%NH4_Urea_m2(3)=NH4_Urea/3;




%========================================================================== 
% Used for testing
% % Delete
% IN_Nit_Fertilizer(day)=4.5;
% IN_Amm_Fertilizer(day)=4.5;
% IN_Organic_Fertilizer(day)=9;
% PPT_mm(day)=1;
% Urea=0;
% Cb=500;
% Ts(1)=15;
% sm(1)=0.3783;
% smp(1)=18./1000;
% % Delete until here
% 
% 




