function [kl, kd, kh] = CN_Compute_kl_kd_kh (PARAMS, VARIABLES, FORCING...
    , VERTSTRUC, SWITCHES, CONSTANTS, vwc_layer_day_point, ice_layer_day_point, temp_layer_day_point...
    , TBMLa, rootfr, WeithedAveCNabove, WeithedAveCNbelow, WeithedAveCNcrRatio, RootDeathWhenH, soiltype, choose_crop);
% function [kl, kd, kh] = CN_Compute_kl_kd_kh (PARAMS, VARIABLES, FORCING...
%    , VERTSTRUC, SWITCHES, CONSTANTS, vwc_layer_day_point, temp_layer_day_point...
%    , TBMLa, rootfr, CNveg, CNveg_root, soiltype);
 
%=========================================================================
% This code solves the linear system of equations (4) (8) and (11) in 
% Porportato (2003) under steady state conditions.
%
%------------------------- Input Variables -------------------------------
%         ADD           % C added to that layer [gC/m^3/d] 
%         Cb            % Total C concetration in the biomass pool for that layer % [gC/m^3]
%         Cl            % Total C concetration in the litter pool for that layer % [gC/m^3]
%         Ch            % Total C concetration in the humus pool for that layer % [gC/m^3]
%         fds           % Nondimensional factor that describes soil
%                         moisture effects on decomposition
%         rh            % Isohumic coefficient that represents the fraction
%                         of decomposition litter that undergoes humification
%         rr            % Portion of decomposing carbon that is lost by
%                        respiration
%         phi           % Nondimensional factor that accounts for a possible reduction of the decomposition rate when the litte   
%
%------------------------- Output Variables ------------------------------
%         kd            % Value that defines death of microbial biomass
%         kl            % Value that defines the rate of decomposition of
%                       litter
%         kh            % Value that defines the rate of decomposition of
%                       humus
%=========================================================================
% Cl, Cb, Ch, and CNl, and CNh
Cl=VARIABLES.Cl;
Cb=VARIABLES.Cb;
Ch=VARIABLES.Ch;
CNl=VARIABLES.CNl;
CNb=PARAMS.CN.CNb;
CNh=PARAMS.CN.CNh;

% Parameters 
rr=PARAMS.CN.rr;
%CR_ratio=PARAMS.CN.CR_ratio;
CR_ratio=WeithedAveCNcrRatio;
nspecies=PARAMS.CanStruc.nspecies;
nl_soil=PARAMS.nl_soil;
dz_mm=VERTSTRUC.dzsmm; 

% TBMLa
TBMLa=FORCING.TBMla;
TBMLaNoLitterBack=FORCING.TBMLaNoLitterBack;
TBMLaOnlyLitterBack=FORCING.TBMLaOnlyLitterBack;
    
% For corn and soy, without LitterBack!, cr factor will be multiplied later
TBMLa_root=TBMLaNoLitterBack;

% Average the entire variables.
%ave_vwc=mean(vwc_layer_day_point+ice_layer_day_point,2);
ave_vwc=mean(vwc_layer_day_point,2);
ave_temp=mean(temp_layer_day_point,2);
ave_TBMLa=mean(TBMLa);
% For corn and soybean
ave_TBMLa_root=mean(TBMLa_root);
ave_CNveg=WeithedAveCNabove;
ave_CNveg_root=WeithedAveCNbelow;
% For corn and soy
ave_rootdeath=mean(RootDeathWhenH);

sm=ave_vwc;
Ts=ave_temp;
TBMla=ave_TBMLa;
CNvegtation=ave_CNveg;
CNvegtation_root=ave_CNveg_root;

% Average the entire
rhmin=PARAMS.CN.rhmin;
rh = min(rhmin, CNh./CNl);

% Soil layer
nl_soil=PARAMS.nl_soil;

% fTd, Paul, K.(2001)
% for i=1:nl_soil
% %    fTd(i) = exp((3.36*(Ts(i)-40))/(Ts(i)+31.79));
%     fTd(i) = 0.125 * exp(0.07*Ts(i));
% end
% indTP = fTd < 0;
% fTd(indTP) = 0;
% indTP = fTd > 1;
% fTd(indTP) = 1;
% fTd = fTd(:);

indT = Ts > 25;
fTd(indT) = 1;              
fTd(~indT) = 1.9.^((Ts(~indT)-25)/10); 
fTd = fTd(:);

% fSd, Paul, K.(2001)
%  / 1/[1+ 4. exp(-6.RWC)] / water content relative to 
% the upper and lower limit of water observed
%in the field (i.e. RWC=(W-LL)/(UL-LL)).
%for i=1:nl_soil
%    fSd(i) = 1/(1+4*exp(-6*((sm(i)-0)/(0.501-0)))); 
%end
%fSd = fSd(:);

% Calculate soil matric potential from soil moisture
for i=1:nl_soil
   if soiltype(i) == 9203 % sand
       res_sm=0.058; % [cm3/cm3]
%       sat_sm=0.37; % [cm3/cm3]
       sat_sm=0.4370001; % [cm3/cm3]
       alpha=0.035;  % [1/cm]
       n=3.19;
       BulkDens=1.69; % [g/cm3]
       FieldCapa=0.12; % [-]
       ClayContent=0.03;
   end
   if soiltype(i) == 8107 % loamy sand
       res_sm=0.074;
%       sat_sm=0.39;
       sat_sm=0.4370001;
       alpha=0.035;
       n=2.39;
       BulkDens=1.63; % [g/cm3]
       FieldCapa=0.14; % [-]
       ClayContent=0.07;
   end
   if soiltype(i) == 6510 % sandy loam
       res_sm=0.067;
%       sat_sm=0.37;
       sat_sm=0.4530001;
       alpha=0.021;
       n=1.61;
       BulkDens=1.51; % [g/cm3]
       FieldCapa=0.23; % [-]
       ClayContent=0.10;
   end
   if soiltype(i) == 4218 % loam
       res_sm=0.083;
       %sat_sm=0.46;
       sat_sm=0.463;
       alpha=0.025;
       n=1.31;
       BulkDens=1.43; % [g/cm3]
       FieldCapa=0.26; % [-]
       ClayContent=0.18;
   end
   
   if soiltype(i) == 2015 % silt loam
       res_sm=0.061;
%       sat_sm=0.43;
       sat_sm= 0.501;
       alpha=0.012;
       n=1.39;
       BulkDens=1.38; % [g/cm3]
       FieldCapa=0.3; % [-]
       ClayContent=0.15;
   end
   if soiltype(i) == 6027 % sandy clay loam
       res_sm=0.086;
%       sat_sm=0.40;
       sat_sm=0.3980001;
       alpha=0.033;
       n=1.49;
       BulkDens=1.40; % [g/cm3]
       FieldCapa=0.33; % [-]
       ClayContent=0.27;
   end
   if soiltype(i) == 3234 % clay loam
       res_sm=0.129;
%       sat_sm=0.47;
       sat_sm=0.464;
       alpha=0.030;
       n=1.37;
       BulkDens=1.30; % [g/cm3]
       FieldCapa=0.335; % [-]
       ClayContent=0.34;
   end
   if soiltype(i) == 933 % silty clay loam
       res_sm=0.098;
%       sat_sm=0.55;
       sat_sm=0.471;
       alpha=0.027;
       n=1.41;
       BulkDens=1.27; % [g/cm3]
       FieldCapa=0.34; % [-]
       ClayContent=0.33;
   end
   if soiltype(i) == 1045 % silty clay
       res_sm=0.163;
%       sat_sm=0.47;
       sat_sm=0.4790001;
       alpha=0.023;
       n=1.39;
       BulkDens=1.21; % [g/cm3]
       FieldCapa=0.36; % [-]
       ClayContent=0.45;
   end
   if soiltype(i) == 2060 % clay
       res_sm=0.102;
%       sat_sm=0.51;
       sat_sm=0.4750001;
       alpha=0.021;
       n=1.20;
       BulkDens=1.25; % [g/cm3]
       FieldCapa=0.36; % [-]
       ClayContent=0.60;
   end
   
   % For matric_potential
   if sm(i) < res_sm
       real_sm(i)=res_sm+0.000000000000001;;
   end
   matric_potential(i)=(1/alpha)*((1/(((sm(i)-res_sm)/(sat_sm-res_sm))^(1/(1-1/n))))-1)^(1/n); % [cm]
   
   % porosity
   porosity(i)=sat_sm;
   
   FieldCapacity(i)=FieldCapa;
   FieldCapacity=FieldCapacity(:);
end
% http://www.convertunits.com/from/cm+H2O/to/megapascal
smp=matric_potential.*0.0000980665; % [cm] to [MPa]

%smp = smp * 9.8066e-06; % First convert to MPa
fSd = log(7.58./abs(smp))/log(7.58/0.01);
indWP = fSd < 0;
fSd(indWP) = 0;
indWP = fSd > 1;
fSd(indWP) = 1;

fSd = fSd(:);
porosity=porosity(:);
if SWITCHES.CN.PoporatofSdfNd_on == 1
    indS = sm(1:end-1) > FieldCapacity; %susana edit
    fSd(indS) = (FieldCapacity(indS)./porosity(indS))./(sm(indS)./porosity(indS));
    fSd(~indS) = (sm(~indS)./porosity(~indS))./(FieldCapacity(~indS)./porosity(~indS));
    fSd=fSd(:);
    
    indSn = sm(1:end-1) > FieldCapacity;
    fSn(indSn) = (1-(sm(indSn)./porosity(indSn)))./(1-(FieldCapacity(indSn)./porosity(indSn)));
    fSn(~indSn) = (sm(~indSn)./porosity(~indSn))./(FieldCapacity(~indSn)./porosity(~indSn));
    fSn=fSn(:);    
    fSn(fSn<0.000001)=0;
end


% phi assumed to be 1
for i=1:nl_soil
    phi(i)=1;
end
phi=phi(:);

% In order to calculate ADD, call 'CN_bioturbation' with average values
% and 'CN_addlitter' with average values
if SWITCHES.CN.Bioturbation;
    [VARIABLES] = CN_bioturbation_with_Mean (PARAMS, VARIABLES...
        , CONSTANTS, FORCING, VERTSTRUC, Cl, TBMla, CNvegtation, CNvegtation_root);

    % For corn and soybean
    [ADD ADD_net] = CN_addlitter_with_Mean (FORCING, PARAMS, VERTSTRUC...
        , VARIABLES, CONSTANTS, SWITCHES, rootfr, TBMla, CNvegtation, CNvegtation_root, ave_rootdeath, ave_TBMLa_root, CR_ratio);
    %[ADD ADD_net] = CN_addlitter_with_Mean (FORCING, PARAMS, VERTSTRUC...
    %    , VARIABLES, CONSTANTS, SWITCHES, rootfr, TBMla, CNvegtation, CNvegtation_root);    
    
    % ADD for Bioturbation is ADD_net
    ADD = ADD_net;    
    if SWITCHES.CN.Bioturbation_new == 1
        [VARIABLES] = CN_bioturbation_with_Mean_new (PARAMS, VARIABLES...
            , CONSTANTS, FORCING, VERTSTRUC, SWITCHES, Cl, Cb, Ch, TBMla, CNvegtation, CNvegtation_root);
        
        [ADD ADD_net] = CN_addlitter_with_Mean_new (FORCING, PARAMS, VERTSTRUC...
            , VARIABLES, CONSTANTS, SWITCHES, rootfr, TBMla, CNvegtation, CNvegtation_root, ave_rootdeath, ave_TBMLa_root, CR_ratio);
        
        ADD = ADD_net;
        
        
        
    end
else
    % Above ground
    ADDa = TBMla; % Total biomass Drop in this time step as litter from above [gr / m^2]     
    %dk ADDa = ADDa * 86400/dtime; % Convert from [gr/m^2/dtime s] to [gr/m^2/d] 
    ADDa = ADDa;

    % Below ground
    ADDb = ADDa*CR_ratio; % Total biomass Drop in this time step as drop from below [gr / m^2/d]  [gC/m^2/d]
    %CNadd(1) = CNveg(1,timestep);
    %CNadd(2) = CNveg(2,timestep);
    %CNadd(3) = CNveg(3,timestep);
    CNadd(1) = CNvegtation_root;
    CNadd(2) = CNvegtation_root;
    CNadd(3) = CNvegtation_root;
    % Litter from plant below ground is distributed based on root distribution

    for jj=1:1:nspecies
        % add litter from below ground
        ADD_r(:,jj)=ADDb(jj).*rootfr./(dz_mm/1000); 
        % add litter above ground is distributed uniformly in first layers     
        ADD_l(:,jj)=zeros(nl_soil,1); 
        ADD_l(1,jj)= ADDa(jj)/(dz_mm(1)/1000);
        %dk if SWITCHES.CN.Bioturbation
        %dk ADD_l(1,jj)= ADD_l(1,jj) - ADDa(jj)/(dz_mm(1)/1000)*biofactor;
        %dk end       
        %ADD_l(2,jj)= ADDa(jj)/2/(dz_mm(2)/1000);         
        % Supersposition of litter from both
        ADDtotal(:,jj)=ADD_r(:,jj)+ADD_l(:,jj);
    end
    % the demention of ADD is [nl_soil, 1]
    ADD = sum(ADDtotal,2);
    
end

    % Test
%     ADD=zeros(nl_soil,1);
%     ADD(1)=TBMla;
%     ADD=ADD+rootfr.*ave_rootdeath+mean(TBMLaNoLitterBack).*rootfr.*CR_ratio;

% solution of Ak = B 
if SWITCHES.CN_type == 1
    for i=1:nl_soil        
        % Defines matrix A
        A(1,1) = Cb(i); 
        A(1,2) = -phi(i)*fSd(i)*fTd(i)*Cb(i)*Cl(i); 
        A(1,3) = 0;

        A(2,1) = 0;
        A(2,2) = (rh(i)*phi(i)*fSd(i)*fTd(i)*Cb(i)*Cl(i));
        A(2,3) = -(phi(i)*fSd(i)*fTd(i)*Cb(i)*Ch(i));

        A(3,1) = -Cb(i);
        A(3,2) = (1-rh(i)-rr)*(phi(i)*fSd(i)*fTd(i)*Cb(i)*Cl(i));
        A(3,3) = (1-rr)*(phi(i)*fSd(i)*fTd(i)*Cb(i)*Ch(i));

        %Defines matrix B
        B(1,1) =-ADD(i);
        B(2,1) = 0;
        B(3,1) = 0;
        if SWITCHES.CN.Bioturbation_new == 1
            B(2,1) = -(VARIABLES.Ch_bio_in(i)-VARIABLES.Ch_bio_out(i));
            B(3,1) = -(VARIABLES.Cb_bio_in(i)-VARIABLES.Cb_bio_out(i));
        end
        
        %Compute the vector
        x=A^(-1)*B;

        kd_transition(i)=x(1);
        kl_transition(i)=x(2);
        kh_transition(i)=x(3);
    end
        kd=kd_transition(:);
        kl=kl_transition(:);
        kh=kh_transition(:);
        
elseif SWITCHES.CN_type == 0

    for i=1:nl_soil
        % Defines matrix A
        A(1,1) = -phi(i)*fSd(i)*fTd(i)*Cb(i)*Cl(i); 
        A(1,2) = Cb(i);

        A(2,1) = phi(i)*fSd(i)*fTd(i)*Cb(i)*Cl(i)*(1-rr); 
        A(2,2) = -Cb(i);


        %Defines matrix B
        B(1,1) =-ADD(i);
        B(2,1) = 0;

        %Compute the vector
        x=A^(-1)*B;

        kl=x(1);
        kd=x(2);
        kh =nan;

        kl_transition(i)=kl;
        kd_transition(i)=kd;
    end
    
    kl=kl_transition(:);
    kd=kd_transition(:);
    kh=nan(nl_soil,1);
end  