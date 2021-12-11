%This function gets the concentrations from Crunch to send back to REWT

function [pH_crunch,Cb_aq,Acetic_acid,C6H12O6_crunch,...
                NH4_crunch,NO3_crunch,N2_crunch,O2_crunch,...
                CO2_crunch,HCO3_crunch,CO32m_crunch,Hplus_crunch,...
                Cb_crunch,NH3_crunch,... CO2gas_crunch,O2gas_crunch,...
                SiO2_crunch, Al3p_crunch,Mg2p_crunch,Ca2p_crunch,...
                Kp_crunch,Fe2p_crunch,Fe3p_crunch,TiO4H4_crunch,Nap_crunch,...
                Albite_crunch, Anatase_crunch, Calcite_crunch,Chamosite_crunch,...
                Anthophyllite_crunch, Kaolinite_crunch, Microcline_crunch,...
                Muscovite_crunch, Quartz_crunch,Goethite_crunch,...
                Lepidocrocite_crunch,Epidote_crunch, Magnetite_crunch] = backfromcrunch()

fid = fopen('influent1layer.out'); % 
crunchoutput = textscan(fid, '%s', 'Delimiter','\t', 'CollectOutput', true,'HeaderLines',2);
fclose(fid);
splitcells = regexp(crunchoutput{:,:}, '\s+', 'split');
splitcells = vertcat(splitcells{:});
crunchresults = str2double(splitcells);
%find where crunch time step closest to 1 day
diff = abs(1-crunchresults(:,1));
rr = find(diff == min(diff));
crunchresults = crunchresults(rr,:);

% columns: VARIABLES = "Time (days)" , "pH", "H+", "C5H7O2NO2(aq)", "Acetic_acid(aq)", "C6H12O6", "NH4+", "NO3-", "N2(aq)", "O2(aq)", "CO2(aq)", "SiO2(aq)", "Al+++", "Mg++", "Ca++", "K+", "Fe++", "Ti(OH)4(aq)", "Na+", "H+", "C5H7O2NO2(aq)", "Acetic_acid(aq)", "C6H12O6", "NH4+", "NO3-", "N2(aq)", "O2(aq)", "CO2(aq)", "SiO2(aq)", "Al+++", "Mg++", "Ca++", "K+", "Fe++", "Ti(OH)4(aq)", "Na+", "OH-", "NH3(aq)", "CO3--", "HCO3-", "
%                            1          2     3         4                 5                 6         7       8       9        10       11          12          13      14     15     16      17       18           19     20         21                  22            23       24      25       26       27         28          29         30       31      32    33      34       35           36    37 
                          
pH_crunch = crunchresults(1,2);
Hplus_crunch = crunchresults(1,3); %total concentration
Cb_aq = crunchresults(1,4); %total concentration
Acetic_acid = crunchresults(1,5); %total concentration
C6H12O6_crunch = crunchresults(1,6); %total concentration
NH4_crunch = 10.^crunchresults(1,24); % individual concentration
NO3_crunch = crunchresults(1,8); %total concentration
N2_crunch = crunchresults(1,9); %total concentration
O2_crunch = crunchresults(1,10); %total concentration
CO2_crunch = 10.^crunchresults(1,28); % individual concentration
SiO2_crunch = crunchresults(1,12); %total concentration
Al3p_crunch = crunchresults(1,13); %total concentration
Mg2p_crunch = crunchresults(1,14); %total concentration
Ca2p_crunch = crunchresults(1,15); %total concentration
Kp_crunch = crunchresults(1,16); %total concentration
Fe2p_crunch = crunchresults(1,17); %total concentration
Fe3p_crunch = crunchresults(1,18); %total concentration
TiO4H4_crunch = crunchresults(1,19); %total concentration
Nap_crunch = crunchresults(1,20); %total concentration
NH3_crunch = 10.^crunchresults(1,38); % individual concentration
CO32m_crunch = 10.^crunchresults(1,39); % individual concentration
HCO3_crunch = 10.^crunchresults(1,40); % individual concentration


% For Microbial Biomass
fid = fopen('volume2.tec');
minoutput = textscan(fid, '%s', 'Delimiter','', 'CollectOutput', true,'HeaderLines',3);
splitcellsmin = regexp(minoutput{:,:}, '\s+', 'split');
splitcellsmin = vertcat(splitcellsmin{:});
crunchresultsmin = str2double(splitcellsmin);

Cb_crunch = crunchresultsmin(1,4);
Albite_crunch = crunchresultsmin(1,5); 
Anatase_crunch = crunchresultsmin(1,6);
Calcite_crunch = crunchresultsmin(1,7); 
Chamosite_crunch = crunchresultsmin(1,8); 
Anthophyllite_crunch = crunchresultsmin(1,9);
Kaolinite_crunch = crunchresultsmin(1,10);
Microcline_crunch=  crunchresultsmin(1,11); 
Muscovite_crunch = crunchresultsmin(1,12); 
Quartz_crunch = crunchresultsmin(1,13); 
Goethite_crunch = crunchresultsmin(1,14);
Lepidocrocite_crunch= crunchresultsmin(1,15); 
Epidote_crunch = crunchresultsmin(1,16); 
Magnetite_crunch = crunchresultsmin(1,17); 

%For gases
% fid = fopen('gases2.out'); % 
% crunchoutput = textscan(fid, '%s', 'Delimiter','\t', 'CollectOutput', true,'HeaderLines',2);
% fclose(fid);
% splitcells = regexp(crunchoutput{:,:}, '\s+', 'split');
% gases = splitcells{2,1};
% crunchresultsgas = str2double(gases);
% CO2gas_crunch = crunchresultsgas(1,2);
% O2gas_crunch = crunchresultsgas(1,3);
end