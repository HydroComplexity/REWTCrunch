function [LCH_solute_Top_Up_m3,LCH_solute_Bottom_Up_m3,LCH_solute,LCH_solute_m2] = fixnegative_solute(solute,LCH_solute_Top_Up_m3,LCH_solute_Bottom_Down_m3,dz_mm, nl_soil)
%==========================================================================  
% Solve negative 
% All of the variables are estimated based on previous time step. This
% means the sum of subtraction variables should not be greater than
% previous state
%
% Do not touch Organic matter related variables
%
% Originally written by Dongkook Woo, UIUC, 2014
% Modified by Susana Roque-Malo, UIUC, 2020
% All rights reserved!
%----------------------------------------------------------------

% 1. Change units to /m2
% 2. Sum all concentrations leaving layer
% 3. compare with concentration in previous time step
% 4. if #4 negative, redistribute leaching to fix
%---------------------------------------------------------------------

%calculation error: omega
omega = 10^-14;

% All variable should be in /m2
solute_m2=solute.*(dz_mm./1000); 

LCH_solute_Top_Up_m2 = LCH_solute_Top_Up_m3.*(dz_mm./1000);
LCH_solute_Bottom_Down_m2 = LCH_solute_Bottom_Down_m3.*(dz_mm./1000);

Check_solute=zeros(nl_soil,1);
for i=1:nl_soil
    % Sum of the subtraction variables
    TotalLoss(i) = LCH_solute_Top_Up_m2(i) + LCH_solute_Bottom_Down_m2(i);
    if TotalLoss(i) ~= 0
        % Check the negativeness
        Check_solute(i)=solute_m2(i)-TotalLoss(i); 
     
        % Redistribute the values
        if Check_solute(i) < 0
            LCH_solute_Top_Up_m2(i)=LCH_solute_Top_Up_m2(i)+Check_solute(i)*(LCH_solute_Top_Up_m2(i)/TotalLoss(i))-(10^(-20));
            if abs(LCH_solute_Top_Up_m2(i)) < omega
                LCH_solute_Top_Up_m2(i) = 0;
            end
            LCH_solute_Bottom_Down_m2(i)=LCH_solute_Bottom_Down_m2(i)+Check_solute(i)*(LCH_solute_Bottom_Down_m2(i)/TotalLoss(i))-(10^(-20));
            if abs(LCH_solute_Top_Up_m2(i)) < omega
                LCH_solute_Top_Up_m2(i) = 0;
            end
        end
    end
end

% Write back
LCH_solute_Top_Up_m3=LCH_solute_Top_Up_m2./(dz_mm./1000);
LCH_solute_Bottom_Up_m2=zeros(nl_soil,1);
LCH_solute_Bottom_Up_m2(1:nl_soil-1)=-LCH_solute_Top_Up_m2(2:nl_soil);
LCH_solute_Bottom_Up_m3=LCH_solute_Bottom_Up_m2./(dz_mm./1000);

LCH_solute_Bottom_Down_m3=LCH_solute_Bottom_Down_m2./(dz_mm./1000);
LCH_solute_Top_Down_m2=zeros(nl_soil,1);
LCH_solute_Top_Down_m2(2:nl_soil)=-LCH_solute_Bottom_Down_m2(1:nl_soil-1);
LCH_solute_Top_Down_m3=LCH_solute_Top_Down_m2./(dz_mm./1000);

% Total Leaching 
LCH_solute=LCH_solute_Top_Up_m3+LCH_solute_Bottom_Up_m3+LCH_solute_Bottom_Down_m3+LCH_solute_Top_Down_m3;
LCH_solute_m2=LCH_solute.*(dz_mm./1000);

end