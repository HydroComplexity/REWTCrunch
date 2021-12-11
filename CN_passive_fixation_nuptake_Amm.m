function [UP_amm_active, Fix_Amm] = CN_passive_fixation_nuptake_Amm (nl_soil, Amm_DEM, SUP_amm, Amm_ku, Amm, NFr);
UP_amm_active = NaN(nl_soil,1); %SUSANA EDIT
Fix_Amm = NaN(nl_soil,1); %SUSANA EDIT

% for Amm
for i=1:nl_soil    
    % for Amm
    Amm_CaseEqn31(i)=Amm_DEM(i)-SUP_amm(i);
    if Amm_CaseEqn31(i) < 0 || Amm_CaseEqn31(i) == 0
        UP_amm_active(i) = 0;  % [g/m2/d]
        Fix_Amm(i) = 0;  % [g/m2/d]
    end
    if Amm_ku(i) < Amm_CaseEqn31(i)
        rup=(Amm_DEM(i)-SUP_amm(i)-SUP_amm(i)*NFr)/(Amm_ku(i)*Amm(i)*(1+NFr));
        if rup < 0
            rup=0;
        end
        if rup > 1
            rup=1;
        end
        if rup == nan
            rup=0;
        end
        
        UP_amm_active(i) = (Amm_ku(i)*Amm(i))*rup;
        % Ensure UP < N
        %IndNeg=(UP_amm_active(i)-Amm_track_m2(i))>0;
        %if IndNeg > 0
        %    UP_amm_active(i)=Amm_track_m2(i);
        %end
        
        %Fix_Amm(i) = min((Amm_ku(i)*Amm(i)*rup+SUP_amm(i))*NFr,Amm_DEM(i)-SUP_amm(i)-Amm_ku(i)*Amm(i)*rup);
        Fix_Amm(i) = min((UP_amm_active(i)+SUP_amm(i))*NFr,Amm_DEM(i)-SUP_amm(i)-UP_amm_active(i));
    end
    if Amm_ku(i) > Amm_CaseEqn31(i) && Amm_CaseEqn31(i) >0
        rup=(Amm_DEM(i)-SUP_amm(i)-SUP_amm(i)*NFr)/((Amm_DEM(i)-SUP_amm(i))*(1+NFr));
        if rup < 0
            rup=0;
        end
        if rup > 1
            rup=1;
        end
        if rup == nan
            rup=0;
        end
        
        UP_amm_active(i) = (Amm_DEM(i)-SUP_amm(i))*rup;
        % Ensure UP < N
        %IndNeg=(UP_amm_active(i)-Amm_track_m2(i))>0;
        %if IndNeg > 0
        %    UP_amm_active(i)=Amm_track_m2(i);
        %end
        
        %Fix_Amm(i) = min(((Amm_DEM(i)-SUP_amm(i))*rup+SUP_amm(i))*NFr,Amm_DEM(i)-SUP_amm(i)-(Amm_DEM(i)-SUP_amm(i))*rup);
        Fix_Amm(i) = min((UP_amm_active(i)+SUP_amm(i))*NFr,Amm_DEM(i)-SUP_amm(i)-UP_amm_active(i));
    end
end