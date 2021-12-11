function [UP_nit_active, Fix_Nit] = CN_passive_fixation_nuptake_Nit (nl_soil, Nit_DEM, SUP_nit, Nit_ku, Nit, NFr);
UP_nit_active = NaN(nl_soil,1); %SUSANA EDIT
Fix_Nit = NaN(nl_soil,1); %SUSANA EDIT
% for Nit
for i=1:nl_soil
    Nit_CaseEqn31(i)=Nit_DEM(i)-SUP_nit(i);
    if Nit_CaseEqn31(i) < 0 || Nit_CaseEqn31(i) == 0
        UP_nit_active(i) = 0;  % [g/m2/d]
        Fix_Nit(i) = 0;  % [g/m2/d]
    end
    if Nit_ku(i) < Nit_CaseEqn31(i)
        rup=(Nit_DEM(i)-SUP_nit(i)-SUP_nit(i)*NFr)/(Nit_ku(i)*Nit(i)*(1+NFr));
        if rup < 0
            rup=0;
        end
        if rup > 1
            rup=1;
        end
        if rup == nan
            rup=0;
        end
        
        UP_nit_active(i) = (Nit_ku(i)*Nit(i))*rup;
        % Ensure UP < N
        %IndNeg=(UP_nit_active(i)-Nit_track_m2(i))>0;
        %if IndNeg > 0
        %    UP_nit_active(i)=Nit_track_m2(i);
        %end
        
        %Fix_Nit(i) = min((Nit_ku(i)*Nit(i)*rup+SUP_nit(i))*NFr,Nit_DEM(i)-SUP_nit(i)-Nit_ku(i)*Nit(i)*rup);
        Fix_Nit(i) = min((UP_nit_active(i)+SUP_nit(i))*NFr,Nit_DEM(i)-SUP_nit(i)-UP_nit_active(i));
    end
    if Nit_ku(i) > Nit_CaseEqn31(i) && Nit_CaseEqn31(i) >0
        rup=(Nit_DEM(i)-SUP_nit(i)-SUP_nit(i)*NFr)/((Nit_DEM(i)-SUP_nit(i))*(1+NFr));
        if rup < 0
            rup=0;
        end
        if rup > 1
            rup=1;
        end
        if rup == nan
            rup=0;
        end
        
        UP_nit_active(i) = (Nit_DEM(i)-SUP_nit(i))*rup;
        % Ensure UP < N
        %IndNeg=(UP_nit_active(i)-Nit_track_m2(i))>0;
        %if IndNeg > 0
        %    UP_nit_active(i)=Nit_track_m2(i);
        %end
        
        %Fix_Nit(i) = min(((Nit_DEM(i)-SUP_nit(i))*rup+SUP_nit(i))*NFr,Nit_DEM(i)-SUP_nit(i)-(Nit_DEM(i)-SUP_nit(i))*rup);
        Fix_Nit(i) = min((UP_nit_active(i)+SUP_nit(i))*NFr,Nit_DEM(i)-SUP_nit(i)-UP_nit_active(i));
    end
end