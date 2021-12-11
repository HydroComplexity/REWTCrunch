function [F_D] = CN_denitrification_factor(VARIABLES, VERTSTRUC, smp, Ts, Cb) 

%==========================================================================  
% GET TD, WD, MD to change rate of denitrification in Crunch
% Based on McGill et al., 1981: Phoenix, a model of the dynamics of carbon
% and nitrogen in grassland soils
% Parton et al., 2001: Generalied model for NOx and N2O emissions from
% soils
%
% Written by Dongkook Woo, UIUC, 2014 - SUSANA ROQUE-MALO  Sept 2020
% All rights reserved!

%------------------------- Input Variables --------------------------------
%       smp             % [MPa] Water potential (http://www.convertunits.com/from/cm+H2O/to/megapascal)
%       Ts              % ['C] Soil temperature
%       Cb              % [gC/m3] Carbon concentration in biomass pool 
%       dz_mm           % [mm] Grid depth 
%       sm              % [m3/m3] Soil moisture (total water in soil: fluid water+soild water)
dz_mm=VERTSTRUC.dzsmm;

%------------------------- Parameters -------------------------------------
smp2bar=10;             % Unit convertion, [MPa] to [bar] 

%------------------------- Output Variables -------------------------------
%       TD            % effect of temperature on denitrification   
%       WD            % effect of moisture on denitrification   
%       MD            % effect of microbial biomass on denitrification   
%==========================================================================    

% TD: Effects of temperature on denitricition [-]
TD=zeros(size(Ts,1),1);

for i=1:size(Ts,1)
    if Ts(i) <= 10
        TD(i)=(0.04/14)*Ts(i)+(0.04/14)*4;
        if TD(i)<=0
            TD(i)=0;
        end
        
    elseif (Ts(i) > 10 && Ts(i) <= 15)
        TD(i)=(0.09/5)*Ts(i)-(0.09/5)*10+0.04;
        
    elseif (Ts(i) > 15 && Ts(i) <= 20)
        TD(i)=(0.32/5)*Ts(i)-(0.32/5)*15+0.13;
        
    elseif (Ts(i) > 20 && Ts(i) <= 30)
        TD(i)=(0.22/10)*Ts(i)-(0.22/10)*20+0.45;
        
    elseif (Ts(i) > 30)
        TD(i)=(0.33/10)*Ts(i)-(0.33/10)*30+0.67;
        if TD(i)>=1
            TD(i)=1;
        end
    end
end

% WD: Effects of moisture on denitrification [-]
sbar=abs(smp).*smp2bar; 
WD=(-1/0.7).*sbar+(1/0.7);
WD(WD>1)=1;
WD(WD<0)=0;
WD=WD';

% MD: Biomass of denitrifiers [gC/m2]
Cb_m2=Cb.*(dz_mm./1000);
MD=0.3.*Cb_m2;

F_D = TD.*WD.*MD;


