function [Vol_m2] = CN_volatilization(NH3_m2, Ts, WFPS, Clay);

%==========================================================================  
% Fertilizer application
% - Volatilization
% Li et al., 1992: A model of nitrous oxide evolution from soil driven by
% rainfall events: 1. Model structure and sensitivity
% % Li, 2000: Modeling trace gas emissions from agricultural ecosystems
%
% Written by Dongkook Woo, UIUC, 2014
% All rights reserved!
%
%------------------------- Input Variables --------------------------------
% Ts                    % ['C] Soil temperature
% NH3                   % [gN/m3] Ammonia
% WFPS                  % [-] Water Filled Porosity 
% Clay                  % [-] Clay amount

%------------------------- Parameters -------------------------------------
T0=45;                   % ['C] Reference temperature

%------------------------- Output Variables -------------------------------
% Vol                % [gN/m2/d] Volatilization

%==========================================================================    
% NH3-N in gas
if Ts(1)<0
    Ts(1)=0;
end
NH3_g=NH3_m2*((Ts(1)/T0)^2); % % NH3-N in gas [gN/m2/d]

% Bulk transfer coefficient
Kb=(1-WFPS(1))*(1-Clay(1));

% Volatilization
Vol_m2(1)=Kb*NH3_g; % [gN/m2/d]