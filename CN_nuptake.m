function [SUP_amm, SUP_nit, SUP_amm_all, SUP_nit_all, SUP_amm_active, SUP_nit_active, SUP_amm_passive, SUP_nit_passive,...
    Amm_DEM, Nit_DEM, Fix_Nit, Fix_Amm, DEM_fraction_l_s_g_b, VARIABLES] = ...
    CN_nuptake (VARIABLES, PARAMS, SWITCHES, VERTSTRUC,...
    a_Nit, a_Amm, sm, layeruptake_all, CARBONS, CNratio_leaf, CNratio_stem, CNratio_grain, CNratio_root, rootfr, choose_crop,...
    scenarios, totlail_layer_day_point, Add_Carbon, how_many_year,...
    detect_first_day_of_year, day, iscorn,...
    DEM_fraction_l_s_g_b_rate,Check_year,~,Tot_UP_Nit_m2_print,Amm_cr,Nit_cr)


        nl_soil=PARAMS.nl_soil;%       nl_soil = # soil layers    
        nspecies = PARAMS.CanStruc.nspecies; % number of species
        Nit = VARIABLES.Nit; 
        Amm = VARIABLES.Amm; 
        Amm_mo = VARIABLES.Amm_mo; 
        Amm = Amm_cr;
        Amm_mo = Amm_cr.*a_Amm;
        Nit = Nit_cr;
%        Nit_track_m2_1st=Nit_track.*(VERTSTRUC.dzsmm./1000);        
%        Amm_track_m2_1st=Amm_track.*(VERTSTRUC.dzsmm./1000);
        
        UP_amm_n=zeros(nl_soil,1); %0
        UP_nit_n=zeros(nl_soil,1); %0
        UP_amm_all=zeros(nl_soil,nspecies); %0
        UP_nit_all=zeros(nl_soil,nspecies); %0
        Recycling = SWITCHES.Recycling; % Switch that decides if recycling of nitrate is on or off
        if Recycling
            factorrec = PARAMS.CN.factorrec; % fraction of Nitrate that is recycled
        end
        dz_in_m=VERTSTRUC.dzsmm./1000;
        Fix_Nit = NaN(nl_soil,1); %SUSANA EDIT
        Fix_Amm = NaN(nl_soil,1); %SUSANA EDIT
        
        nfix_rate=PARAMS.nfix_rate;
        NFr=nfix_rate(day);
        
        Nremobilization_rate=PARAMS.Nremobilization_rate;
        Nremobil=SWITCHES.Nremobil;
        Nremobil = 0; %SUSANA OVERWRITE
%*************************************************************************
%       Nitrogen UP take from soil
%*************************************************************************
        % PLANT UPTAKE from soil
        % AMMONIUM
        
        for i=1:1:nspecies
            layeruptn=layeruptake_all(1:end-1,i); % plant water uptake %Susana edit for indeces
            wtind=layeruptn>0;
            wtnind=layeruptn<0;
            wtzind=layeruptn==0;
            % First for the case in which the flow is from the soil to the root
            %UP_amm_n(wtind) = (a_Amm./sm(wtind).*Amm(wtind)).*layeruptn(wtind);% [gr/m^2/d];
            UP_amm_n(wtind) = (Amm_mo(wtind)./sm(wtind)).*layeruptn(wtind);% [gr/m^2/d];
            % When the flow is from the root to the soil. It is assumed is
            % completely mixed inside the root system.
            if Recycling
                waterrec = sum(UP_amm_n(wtind))*factorrec;
                UP_amm_n(wtnind) = - waterrec.*layeruptn(wtnind)/(sum(layeruptn(wtnind)));% [gr/m^2/d]
            else
                UP_amm_n(wtnind) = 0;
            end
            UP_amm_n(wtzind) = 0;

            UP_amm_n=UP_amm_n(:);
            UP_amm_all(:,i) = UP_amm_n;


            % NITRATE
            % First for the case in which the flow is from the soil to the root
            UP_nit_n(wtind) = (a_Nit./sm(wtind).*Nit(wtind)).*layeruptn(wtind);%[gr/m^2/d];
            % When the flow is from the root to the soil. It is assumed is
            % completely mixed inside the root system.
            if Recycling
                waterrec = sum(UP_nit_n(wtind))*factorrec;
                UP_nit_n(wtnind) = - waterrec.*layeruptn(wtnind)/(sum(layeruptn(wtnind)));% [gr/m^2/d]
            else
                UP_nit_n(wtnind) = 0;
            end
            UP_nit_n(wtzind) = 0;

            UP_nit_n=UP_nit_n(:);
            UP_nit_all(:,i) = UP_nit_n;


            clear UP_nit_n;
            clear UP_amm_n;

            UP_nit_active=NaN(size(UP_nit_all));
            UP_amm_active=nan;
            Amm_DEM=nan;
            Nit_DEM=nan;
        end
        SUP_amm = sum(UP_amm_all,2);         % [gr/m^2/d]
        SUP_nit = sum(UP_nit_all,2);         % [gr/m^2/d]
        
        % Ensure UP<N
        %IndNeg=(SUP_amm-Amm_track_m2_1st)>0;
        %SUP_amm(IndNeg)=Amm_track_m2_1st(IndNeg);
        %Amm_track_m2_2nd=Amm_track_m2_1st-SUP_amm;
        %IndNeg=(SUP_nit-Nit_track_m2_1st)>0;
        %SUP_nit(IndNeg)=Nit_track_m2_1st(IndNeg);        
        %Nit_track_m2_2nd=Nit_track_m2_1st-SUP_nit;

        
        
        Active_Nuptake=SWITCHES.Active_Nuptake;
        if Active_Nuptake
            % Nit, and Amm in the unit of [g/m3]
            % carbon content's units from PALMS [kg/m2]
            grain_carbon=CARBONS.grain_carbon_day_point(day)*1000;     % [g/m2]
            leaf_carbon=CARBONS.leaf_carbon_day_point(day)*1000;       % [g/m2]
            rhizome_carbon=CARBONS.rhizome_carbon_day_point(day)*1000; % [g/m2]
            root_carbon=CARBONS.root_carbon_day_point(day)*1000;       % [g/m2]
            stem_carbon=CARBONS.stem_carbon_day_point(day)*1000;       % [g/m2]

            %CNratio=ALL_CNveg(day);
            %CNratio_root=ALL_CNveg_root(day);
            CNleaf=CNratio_leaf(day);
            CNstem=CNratio_stem(day);
            CNgrain=CNratio_grain(day);
            CNroot=CNratio_root(day);

            if day>1
                pre_grain_carbon=CARBONS.grain_carbon_day_point(day-1)*1000;
                pre_leaf_carbon=CARBONS.leaf_carbon_day_point(day-1)*1000;
                pre_rhizome_carbon=CARBONS.rhizome_carbon_day_point(day-1)*1000;
                pre_root_carbon=CARBONS.root_carbon_day_point(day-1)*1000;
                pre_stem_carbon=CARBONS.stem_carbon_day_point(day-1)*1000;

                %pre_CNratio=ALL_CNveg(day-1);
                %pre_CNratio_root=ALL_CNveg_root(day-1);
                pre_CNleaf=CNratio_leaf(day-1);
                pre_CNstem=CNratio_stem(day-1);
                pre_CNgrain=CNratio_grain(day-1);
                pre_CNroot=CNratio_root(day-1);
            else
                pre_grain_carbon=CARBONS.grain_carbon_day_point(day)*1000;
                pre_leaf_carbon=CARBONS.leaf_carbon_day_point(day)*1000;
                pre_rhizome_carbon=CARBONS.rhizome_carbon_day_point(day)*1000;
                pre_root_carbon=CARBONS.root_carbon_day_point(day)*1000;
                pre_stem_carbon=CARBONS.stem_carbon_day_point(day)*1000;

                %pre_CNratio=ALL_CNveg(day);
                %pre_CNratio_root=ALL_CNveg_root(day);
                pre_CNleaf=CNratio_leaf(day);
                pre_CNstem=CNratio_stem(day);
                pre_CNgrain=CNratio_grain(day);
                pre_CNroot=CNratio_root(day);
            end

            F=PARAMS.F;
            d=PARAMS.d;
            n=PARAMS.porosity;
            n=n(:);
            % Calculate ku: dependence of the diffusion process on
            % the soil moisture level [-/d]
            Nit_ku=(a_Nit./(sm))*F.*((sm./n).^d); % [m/d]
            %Amm_ku=(a_Amm./(sm))*F.*((sm./n).^d); % [m/d]
            Amm_ku=(1./(sm))*F.*((sm./n).^d); % [m/d]

            % Calculate DEM based on carbon content change.
            %Above_Current_total_C=grain_carbon+leaf_carbon+stem_carbon; % [g/m2]
            Below_Current_total_C=rhizome_carbon+root_carbon; % [g/m2]

            %Above_Previous_total_C=pre_grain_carbon+pre_leaf_carbon+pre_stem_carbon; % [g/m2]
            Below_Previous_total_C=pre_rhizome_carbon+pre_root_carbon; % [g/m2]

            %Above_C_change=Above_Current_total_C-Above_Previous_total_C;
            leaf_C_change=leaf_carbon-pre_leaf_carbon;
            stem_C_change=stem_carbon-pre_stem_carbon;
            grain_C_change=grain_carbon-pre_grain_carbon;
            Below_C_change=Below_Current_total_C-Below_Previous_total_C;

            % Calibaration
            %if Above_C_change > 0
            %    Above_C_change=(Above_Current_total_C-Above_Previous_total_C)*Add_Carbon; % [g/m2]
            %end
            %if Below_C_change > 0
            %    Below_C_change=(Below_Current_total_C-Below_Previous_total_C)*Add_Carbon; % [g/m2]
            %end

            % in order to have 0 carbon change for every first day of year
            for i=1:how_many_year
                if day == detect_first_day_of_year(i)
                    %Above_C_change=0;
                    leaf_C_change=0;
                    stem_C_change=0;
                    grain_C_change=0;
                    Below_C_change=0;
                end
            end

            % average CN ratio
            %ave_CNratio=(pre_CNratio+CNratio)/2;
            ave_leaf_CNratio=(pre_CNleaf+CNleaf)/2;
            ave_stem_CNratio=(pre_CNstem+CNstem)/2;
            ave_grain_CNratio=(pre_CNgrain+CNgrain)/2;
            ave_CNratio_root=(pre_CNroot+CNroot)/2;

            % Calculate Demand of N
            %if Above_C_change > 0
            %    Above_DEM=(Above_C_change/ave_CNratio).*rootfr; % [g/m2]
            %else
            %    Above_DEM=zeros(nl_soil,1);
            %end
            if leaf_C_change > 0
                leaf_DEM=(leaf_C_change/ave_leaf_CNratio).*rootfr; % [g/m2]
            else
                leaf_DEM=zeros(nl_soil,1);
            end
            if stem_C_change > 0
                stem_DEM=(stem_C_change/(ave_stem_CNratio)).*rootfr; % [g/m2]
            else
                stem_DEM=zeros(nl_soil,1);
            end
            if grain_C_change > 0
                grain_DEM=(grain_C_change/ave_grain_CNratio).*rootfr; % [g/m2]
            else
                grain_DEM=zeros(nl_soil,1);
            end

            if Below_C_change > 0
                Below_DEM=(Below_C_change/ave_CNratio_root).*rootfr; % [g/m2]
            else
                Below_DEM=zeros(nl_soil,1);
            end
            %DEM=Above_DEM+Below_DEM;
            DEM=leaf_DEM+stem_DEM+grain_DEM+Below_DEM;
            
            % For GUI
            VARIABLES.DEM=sum(DEM); % [g/m2]
            
            % DEM fraction
            if sum(DEM) == 0
                DEM_fraction_l_s_g_b = [0 0 0 0];
            else
                DEM_fraction_l_s_g_b=[sum(leaf_DEM) sum(stem_DEM) sum(grain_DEM) sum(Below_DEM)]./sum(DEM);
            end
            
            % Due to low DEM
%             if choose_crop ==4
%                 if day < 365*3
%                     DEM=DEM*3;
%                 end
%                 if day < 365*2
%                     DEM=DEM/3*10;
%                 end
%                 if day < 365
%                     DEM=DEM/3/10*30;
%                 end
%             end

            % Corn and soy
            if choose_crop ==2
                if iscorn == 1 % Corn
                    Nit_DEM=DEM.*0.00;%(0.5/(0.2+0.5)); % The values (0.2 and 0.5) came from D'Odorico et al., 2003
                    Amm_DEM=DEM.*1.00;%(0.2/(0.2+0.5)); % in which, they use constant DEM+-, 0.7 to 0/85 for DEM-
                else % Soy
                    Nit_DEM=DEM.*0.00;%(0.5/(0.2+0.5)); % The values (0.2 and 0.5) came from D'Odorico et al., 2003
                    Amm_DEM=DEM.*1.00;%(0.2/(0.2+0.5)); % in which, they use constant DEM+-, 0.7 to 0/85 for DEM-
                end
            end
            % SG: 0.85, 0.15
            if choose_crop ==3
                Nit_DEM=DEM.*0.9;%(0.5/(0.2+0.5)); % The values (0.2 and 0.5) came from D'Odorico et al., 2003
                Amm_DEM=DEM.*0.1;%(0.2/(0.2+0.5)); % in which, they use constant DEM+-, 0.7 to 0/85 for DEM-
            end
            % MG
            if choose_crop ==4
                Nit_DEM=DEM.*0.45;%0.45; % The values (0.2 and 0.5) came from D'Odorico et al., 2003
                Amm_DEM=DEM.*0.55;%0.55; % in which, they use constant DEM+-, 0.7 to 0/85 for DEM-
            end

            % for NP, I choose Nit_DEM and Amm_DEM based on the paper by Porporato et al. (2003).
            if choose_crop ==1
                %if day>1
                %    pre_totlail=totlail_layer_day_point(day-1);
                %else
                %    pre_totlail=totlail_layer_day_point(day);
                %end
                %totlail=totlail_layer_day_point(day);
                %lai_change=totlail-pre_totlail;

                % 0.8m -> 0.5, 0.2 when rootfr = 0.9760
                % 0.5*(1+1-0.9760)=0.512
                % 0.2*(1+1-0.9760)=0.2048

                %if lai_change > 0
                Nit_DEM=(rootfr.*0.075).*dz_in_m; % 0.5 [g/m3/d] * soil depth or 0.06
                Amm_DEM=(rootfr.*0.000).*dz_in_m; % 0.2 [g/m3/d] * soil depth or 0.01
                %else
                %    Nit_DEM=zeros(nl_soil,1);
                %    Amm_DEM=zeros(nl_soil,1);
                %end
            end

            % Porporato et al. (2003) 's three cases in eqn. 31.
            % UP_amm, UP_nit [gr/m^2/d]
            % Nit, and Amm [g/m3]
            % DEM [g/m2]
            % Nit_ku, Amm_ku [m/d]

            
            [UP_nit_active, Fix_Nit] = CN_passive_fixation_nuptake_Nit (nl_soil, Nit_DEM, SUP_nit, Nit_ku, Nit, NFr);
            [UP_amm_active, Fix_Amm] = CN_passive_fixation_nuptake_Amm (nl_soil, Amm_DEM, SUP_amm, Amm_ku, Amm_mo, NFr);            
            UP_nit_active=UP_nit_active(:);
            UP_amm_active=UP_amm_active(:);
            % Ensure UP<N
            %Nit_track_m2_3rd=Nit_track_m2_2nd-UP_nit_active;
            %Amm_track_m2_3rd=Amm_track_m2_2nd-UP_amm_active;
            
%             % for Nit
%             for i=1:nl_soil
%                 Nit_CaseEqn31(i)=Nit_DEM(i)-SUP_nit(i);
%                 if Nit_CaseEqn31(i) < 0 || Nit_CaseEqn31(i) == 0
%                     UP_nit_active(i) = 0;  % [g/m2/d]
%                     Fix_Nit(i) = 0;  % [g/m2/d]
%                 end
%                 if Nit_ku(i) < Nit_CaseEqn31(i)
%                     rup=(Nit_DEM(i)-SUP_nit(i)-SUP_nit(i)*NFr)/(Nit_ku(i)*Nit(i)*(1+NFr));
%                     if rup < 0
%                         rup=0;
%                     end
%                     if rup > 1
%                         rup=1;
%                     end
%                     if rup == nan
%                         rup=0;
%                     end
% 
%                     UP_nit_active(i) = (Nit_ku(i)*Nit(i))*rup;
%                     Fix_Nit(i) = min((Nit_ku(i)*Nit(i)*rup+SUP_nit(i))*NFr,Nit_DEM(i)-SUP_nit(i)-Nit_ku(i)*Nit(i)*rup);
%                 end
%                 if Nit_ku(i) > Nit_CaseEqn31(i) && Nit_CaseEqn31(i) >0
%                     rup=(Nit_DEM(i)-SUP_nit(i)-SUP_nit(i)*NFr)/((Nit_DEM(i)-SUP_nit(i))*(1+NFr));
%                     if rup < 0
%                         rup=0;
%                     end
%                     if rup > 1
%                         rup=1;
%                     end
%                     if rup == nan
%                         rup=0;
%                     end
% 
%                     UP_nit_active(i) = (Nit_DEM(i)-SUP_nit(i))*rup;
%                     Fix_Nit(i) = min(((Nit_DEM(i)-SUP_nit(i))*rup+SUP_nit(i))*NFr,Nit_DEM(i)-SUP_nit(i)-(Nit_DEM(i)-SUP_nit(i))*rup);
%                 end
% 
%                 % for Amm
%                 Amm_CaseEqn31(i)=Amm_DEM(i)-SUP_amm(i);
%                 if Amm_CaseEqn31(i) < 0 || Amm_CaseEqn31(i) == 0
%                     UP_amm_active(i) = 0;  % [g/m2/d]
%                     Fix_Amm(i) = 0;  % [g/m2/d]
%                 end
%                 if Amm_ku(i) < Amm_CaseEqn31(i)
%                     rup=(Amm_DEM(i)-SUP_amm(i)-SUP_amm(i)*NFr)/(Amm_ku(i)*Amm(i)*(1+NFr));
%                     if rup < 0
%                         rup=0;
%                     end
%                     if rup > 1
%                         rup=1;
%                     end
%                     if rup == nan
%                         rup=0;
%                     end
% 
%                     UP_amm_active(i) = (Amm_ku(i)*Amm(i))*rup;
%                     Fix_Amm(i) = min((Amm_ku(i)*Amm(i)*rup+SUP_amm(i))*NFr,Amm_DEM(i)-SUP_amm(i)-Amm_ku(i)*Amm(i)*rup);
%                 end
%                 if Amm_ku(i) > Amm_CaseEqn31(i) && Amm_CaseEqn31(i) >0
%                     rup=(Amm_DEM(i)-SUP_amm(i)-SUP_amm(i)*NFr)/((Amm_DEM(i)-SUP_amm(i))*(1+NFr));
%                     if rup < 0
%                         rup=0;
%                     end
%                     if rup > 1
%                         rup=1;
%                     end
%                     if rup == nan
%                         rup=0;
%                     end
% 
%                     UP_amm_active(i) = (Amm_DEM(i)-SUP_amm(i))*rup;
%                     Fix_Amm(i) = min(((Amm_DEM(i)-SUP_amm(i))*rup+SUP_amm(i))*NFr,Amm_DEM(i)-SUP_amm(i)-(Amm_DEM(i)-SUP_amm(i))*rup);
%                 end
%             end
            SUP_nit_active=UP_nit_active(:); % [g/m2/d]
            SUP_amm_active=UP_amm_active(:); % [g/m2/d]

            SUP_nit_passive=SUP_nit(:); % [g/m2/d]
            SUP_amm_passive=SUP_amm(:); % [g/m2/d]

            % N Uptake from soil = passive upatke + active uptake
            SUP_nit=(SUP_nit_passive+SUP_nit_active); % [g/m2/d]
            SUP_amm=(SUP_amm_passive+SUP_amm_active); % [g/m2/d]

            SUP_nit_all=SUP_nit; % [g/m2/d]
            SUP_amm_all=SUP_amm; % [g/m2/d]

            % N Uptake from N fixation
            Fix_Nit=Fix_Nit(:);
            Fix_Amm=Fix_Amm(:);
            
        end
        
        
        if Nremobil
            for i=1:how_many_year
                if day == detect_first_day_of_year(i)
                    if Check_year == 1
                        TUPr_yrAmm_pool=0;
                        TUPr_yrNit_pool=0;
                        
                        VARIABLES.TUPr_yrAmm_pool=TUPr_yrAmm_pool;
                        VARIABLES.TUPr_yrNit_pool=TUPr_yrNit_pool;
                    else
                        DEM_fraction_Abov_p_below=(1-DEM_fraction_l_s_g_b_rate(:,4))';
                        TUPr_yrAmm_pool=sum(sum(Tot_UP_Amm_m2_print(:,detect_first_day_of_year(Check_year-1):detect_first_day_of_year(Check_year)-1)).* ...
                            (DEM_fraction_Abov_p_below(:,detect_first_day_of_year(Check_year-1):detect_first_day_of_year(Check_year)-1)).* ...
                            (Nremobilization_rate(:,detect_first_day_of_year(Check_year-1):detect_first_day_of_year(Check_year)-1))); % g/m2
                        TUPr_yrNit_pool=sum(sum(Tot_UP_Nit_m2_print(:,detect_first_day_of_year(Check_year-1):detect_first_day_of_year(Check_year)-1)).* ...
                            (DEM_fraction_Abov_p_below(:,detect_first_day_of_year(Check_year-1):detect_first_day_of_year(Check_year)-1)).* ...
                            (Nremobilization_rate(:,detect_first_day_of_year(Check_year-1):detect_first_day_of_year(Check_year)-1)));
 
                        VARIABLES.TUPr_yrAmm_pool=TUPr_yrAmm_pool;
                        VARIABLES.TUPr_yrNit_pool=TUPr_yrNit_pool;
                    end
                end
            end
            
            
            TUPr_yrAmm_pool=VARIABLES.TUPr_yrAmm_pool;
            TUPr_yrNit_pool=VARIABLES.TUPr_yrNit_pool;
            UPr_amm=zeros(nl_soil,1);
            UPr_nit=zeros(nl_soil,1);
            
            % Ammonium
            if TUPr_yrAmm_pool > 0
                if TUPr_yrAmm_pool >= sum(Amm_DEM)
                    UPr_amm=Amm_DEM;
                    TUPr_yrAmm_pool=TUPr_yrAmm_pool-sum(Amm_DEM); % g/m2
                    
                    % Update active, passive, fixation uptakes
                    SUP_amm_active=zeros(nl_soil,1); % [g/m2/d]
                    SUP_amm_passive=zeros(nl_soil,1); % [g/m2/d]
                    
                    % N Uptake from soil = passive upatke + active uptake
                    SUP_amm=(SUP_amm_passive+SUP_amm_active); % [g/m2/d]
                    SUP_amm_all=SUP_amm; % [g/m2/d]
                    
                    % N Uptake from N fixation
                    Fix_Amm=zeros(nl_soil,1);
                else
                    UPr_amm=TUPr_yrAmm_pool.*rootfr;
                    TUPr_yrAmm_pool=0; % g/m2
                    
                    % Update active, passive, fixation uptakes
                    Amm_DEM_updated=(sum(Amm_DEM)- sum(UPr_amm)).*rootfr;
                    [UP_amm_active, Fix_Amm] = CN_passive_fixation_nuptake_Amm (nl_soil, Amm_DEM_updated, SUP_amm, Amm_ku, Amm_mo, NFr);            
            
                    SUP_amm_active=UP_amm_active(:); % [g/m2/d]
                    SUP_amm_passive=SUP_amm_passive; % [g/m2/d]
                    
                    % N Uptake from soil = passive upatke + active uptake
                    SUP_amm=(SUP_amm_passive+SUP_amm_active); % [g/m2/d]
                    SUP_amm_all=SUP_amm; % [g/m2/d]
                    
                    % N Uptake from N fixation
                    Fix_Amm=Fix_Amm(:);
                end
            end
            VARIABLES.TUPr_yrAmm_pool=TUPr_yrAmm_pool;
            
            % Nitrate
            if TUPr_yrNit_pool > 0
                if TUPr_yrNit_pool >= sum(Nit_DEM)
                    UPr_nit=Nit_DEM;
                    TUPr_yrNit_pool=TUPr_yrNit_pool-sum(Nit_DEM); % g/m2
                    
                    % Update active, passive, fixation uptakes
                    SUP_nit_active=zeros(nl_soil,1); % [g/m2/d]
                    SUP_nit_passive=zeros(nl_soil,1); % [g/m2/d]
                    
                    % N Uptake from soil = passive upatke + active uptake
                    SUP_nit=(SUP_nit_passive+SUP_nit_active); % [g/m2/d]
                    SUP_nit_all=SUP_nit; % [g/m2/d]
                    
                    % N Uptake from N fixation
                    Fix_Nit=zeros(nl_soil,1);
                else
                    UPr_nit=TUPr_yrNit_pool.*rootfr;
                    TUPr_yrNit_pool=0; % g/m2
                     
                    % Update active, passive, fixation uptakes
                    Nit_DEM_updated=(sum(Nit_DEM)- sum(UPr_nit)).*rootfr;
                    [UP_nit_active, Fix_Nit] = CN_passive_fixation_nuptake_Nit (nl_soil, Nit_DEM_updated, SUP_nit, Nit_ku, Nit, NFr);
                                
                    SUP_nit_active=UP_nit_active(:); % [g/m2/d]
                    SUP_nit_passive=SUP_nit_passive; % [g/m2/d]
                    
                    % N Uptake from soil = passive upatke + active uptake
                    SUP_nit=(SUP_nit_passive+SUP_nit_active); % [g/m2/d]
                    SUP_nit_all=SUP_nit; % [g/m2/d]
                    
                    % N Uptake from N fixation
                    Fix_Nit=Fix_Nit(:);
                end
            end
            VARIABLES.TUPr_yrNit_pool=TUPr_yrNit_pool;
            
        end
 
        % Track
        %Nit_track_m2=Nit_track_m2_1st-SUP_nit_passive-SUP_nit_active;
        %Amm_track_m2=Amm_track_m2_1st-SUP_amm_passive-SUP_amm_active;
        %Nit_track=Nit_track_m2./(VERTSTRUC.dzsmm./1000); 
        %Amm_track=Amm_track_m2./(VERTSTRUC.dzsmm./1000);
        
% %*************************************************************************
% %       Nitrogen Fixation
% %*************************************************************************
% % Assution: 1. if there is N fertilizer, the crop meet N demand from soil,
% % but some propotion come from N fixation. Thus soil N goes back eventhough
% % N uptake remains the same.
% %           2. if there is No N fertilizer, the crop do not meet N demand
% % from soil thus, Just increase N demand by seperating process. Thus, When 
% % you save the whole process, you need to multiply N fixation percentage
% % to the final save point!
% nfix_rate=PARAMS.nfix_rate;
% 
%Fix_Nit=zeros(size(SUP_nit,1),1);
%Fix_Amm=zeros(size(SUP_amm,1),1);
% 
% if SWITCHES.CN.Nfix == 1
%     if choose_crop == 2
%         % soybean case
%         if iscorn == 0
%             for i=1:nl_soil
%                 % Nit
%                 if Nit_DEM(i) <= SUP_nit_active(i) + SUP_nit_passive(i)
%                     Fix_Nit(i)=(SUP_nit_active(i) + SUP_nit_passive(i)).*nfix_rate;
%                     
%                     SSUP_nit_active(i)=(SUP_nit_active(i) + SUP_nit_passive(i) - Fix_Nit(i)).*((SUP_nit_active(i))./(SUP_nit_active(i) + SUP_nit_passive(i)));
%                     SSUP_nit_passive(i)=(SUP_nit_active(i) + SUP_nit_passive(i) - Fix_Nit(i)).*((SUP_nit_passive(i))./(SUP_nit_active(i) + SUP_nit_passive(i)));
%                 end
%                 if Nit_DEM(i) > SUP_nit_active(i) + SUP_nit_passive(i)
%                     Fix_Nit_max(i)=Nit_DEM(i).*nfix_rate;
%                     if SUP_nit_active(i) + SUP_nit_passive(i) <= Nit_DEM(i) - Fix_Nit_max(i) 
%                         Fix_Nit(i)=(SUP_nit_active(i) + SUP_nit_passive(i)).*nfix_rate;
%                         
%                         SSUP_nit_active(i)=SUP_nit_active(i);
%                         SSUP_nit_passive(i)=SUP_nit_active(i);
%                     end
%                     if SUP_nit_active(i) + SUP_nit_passive(i) > Nit_DEM(i) - Fix_Nit_max(i)
%                         Fix_Nit(i)=(SUP_nit_active(i) + SUP_nit_passive(i)).*nfix_rate;
%                         
%                         SSUP_nit_active(i)=(SUP_nit_active(i) + SUP_nit_passive(i) - Fix_Nit(i)).*((SUP_nit_active(i))./(SUP_nit_active(i) + SUP_nit_passive(i)));
%                         SSUP_nit_passive(i)=(SUP_nit_active(i) + SUP_nit_passive(i) - Fix_Nit(i)).*((SUP_nit_passive(i))./(SUP_nit_active(i) + SUP_nit_passive(i)));                
%                     end
%                 end
%                               
%                 % Amm
%                 if Amm_DEM(i) <= SUP_amm_active(i) + SUP_amm_passive(i)
%                     Fix_Amm(i)=(SUP_amm_active(i) + SUP_amm_passive(i)).*nfix_rate;
%                     
%                     SSUP_amm_active(i)=(SUP_amm_active(i) + SUP_amm_passive(i) - Fix_Amm(i)).*((SUP_amm_active(i))./(SUP_amm_active(i) + SUP_amm_passive(i)));
%                     SSUP_amm_passive(i)=(SUP_amm_active(i) + SUP_amm_passive(i) - Fix_Amm(i)).*((SUP_amm_passive(i))./(SUP_amm_active(i) + SUP_amm_passive(i)));
%                 end
%                 if Amm_DEM(i) > SUP_amm_active(i) + SUP_amm_passive(i)
%                     Fix_Amm_max(i)=Amm_DEM(i).*nfix_rate;
%                     if SUP_amm_active(i) + SUP_amm_passive(i) <= Amm_DEM(i) - Fix_Nit_max(i) 
%                         Fix_Nit(i)=(SUP_nit_active(i) + SUP_nit_passive(i)).*nfix_rate;
%                         
%                         SSUP_nit_active(i)=SUP_nit_active(i);
%                         SSUP_nit_passive(i)=SUP_nit_active(i);
%                     end
%                     if SUP_nit_active(i) + SUP_nit_passive(i) > Nit_DEM(i) - Fix_Nit_max(i)
%                         Fix_Nit(i)=(SUP_nit_active(i) + SUP_nit_passive(i)).*nfix_rate;
%                         
%                         SSUP_nit_active(i)=(SUP_nit_active(i) + SUP_nit_passive(i) - Fix_Nit(i)).*((SUP_nit_active(i))./(SUP_nit_active(i) + SUP_nit_passive(i)));
%                         SSUP_nit_passive(i)=(SUP_nit_active(i) + SUP_nit_passive(i) - Fix_Nit(i)).*((SUP_nit_passive(i))./(SUP_nit_active(i) + SUP_nit_passive(i)));                
%                     end
%                 end
%             end
%             %nfix_rate=PARAMS.nfix_rate;
%             
%         end
%     end
% end