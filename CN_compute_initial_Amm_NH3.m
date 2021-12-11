function [Amm, NH3, Amm_mo, Amm_im, VARIABLES] = CN_compute_initial_Amm_NH3(PARAMS, temp_layer_day_point, nl_soil, VARIABLES)

% Caution: Also look at "CN_ammonia_ammonium_exchange"

Ts=mean(temp_layer_day_point,2);
Amm=VARIABLES.Amm;
NH3=VARIABLES.NH3;
pH=VARIABLES.pH;

%%-----------------------------------------------------------------------%%
% Exchange between NH4 and NH3
CN_NH4m_fraction ();

% As can be shown, there are only one N in NH4 and NH3. Thus, we can use
% the fraction of NH4m directly
NH3_NH4=Amm+NH3;

New_NH4=NH3_NH4.*fraction_NH4m; % [gN/m3]
New_NH3=NH3_NH4.*(1-fraction_NH4m); % [gN/m3]

%%-----------------------------------------------------------------------%%
VARIABLES.Amm=New_NH4;
Amm=VARIABLES.Amm;
VARIABLES.NH3=New_NH3;
NH3=VARIABLES.NH3;

% Amm mobile and immobile
VARIABLES.Amm_mo=VARIABLES.Amm.*PARAMS.CN.a_Amm;
Amm_mo=VARIABLES.Amm_mo;
VARIABLES.Amm_im=VARIABLES.Amm.*(1-PARAMS.CN.a_Amm);
Amm_im=VARIABLES.Amm_im;
