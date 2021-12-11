%%-----------------------------------------------------------------------%%
% 1. Load the tiledrainage
% total m3/d in each field, each plot is in 200m x 200m
for i=start_year:end_year
    a='data2NcycleTiledrain';
    b=num2str(i);
    c='_S';
    d=num2str(scenarios);
    e=strcat(a,b,c,d);
    
    f='Daily_sum_drain_NP';
    g='Daily_sum_drain_Corn';
    ff='Daily_sum_drain_MG';
    gg='Daily_sum_drain_SG';
    
    load(e,f);
    Daily_sum_drain_NP=Daily_sum_drain_NP';
    load(e,g);
    Daily_sum_drain_Corn=Daily_sum_drain_Corn';
    load(e,ff);
    Daily_sum_drain_MG=Daily_sum_drain_MG';
    load(e,gg);
    Daily_sum_drain_SG=Daily_sum_drain_SG';
    
    if i == start_year
        % total m3/d in each field, each plot is in 200m x 200m
        Drain_NP=[Daily_sum_drain_NP];
        Drain_Corn=[Daily_sum_drain_Corn];
        Drain_MG=[Daily_sum_drain_MG];
        Drain_SG=[Daily_sum_drain_SG];
    else
        Drain_NP=[Drain_NP;Daily_sum_drain_NP];
        Drain_Corn=[Drain_Corn;Daily_sum_drain_Corn];
        Drain_MG=[Drain_MG;Daily_sum_drain_MG];
        Drain_SG=[Drain_SG;Daily_sum_drain_SG];
    end
    % SUSANA EDIT: NO TILE DRAINAGE.
    Drain_NP = zeros(size(Drain_NP));
    Drain_Corn = zeros(size(Drain_Corn));
    Drain_MG = zeros(size(Drain_MG));
    Drain_SG = zeros(size(Drain_SG));

end

%%-----------------------------------------------------------------------%%
% 2. Load variables from PALMS' outputs
% SUSANA EDIT
if start_year == 2008
    load('noPALMSforcing2008.mat')
    qq_layer_day_point = qq_noPALMS;
    ice_layer_day_point = soilice_noPALMS;
    vwc_layer_day_point = Sm_noPALMS;
    temp_layer_day_point = Ts_noPALMS;
    adupsoil_layer_day_point = adup_noPALMS;
    load('noPALMScorn');
    grain_carbon_day_point = grainC_corn./1000;
    leaf_carbon_day_point = leafC_corn./1000;
    rhizome_carbon_day_point = rhizC_corn./1000;
    root_carbon_day_point = rootC_corn./1000;
    stem_carbon_day_point = stemC_corn./1000;
else
    load('noPALMSforcing2010.mat')
    qq_layer_day_point = qq_noPALMS;
    ice_layer_day_point = soilice_noPALMS;
    vwc_layer_day_point = Sm_noPALMS;
    temp_layer_day_point = Ts_noPALMS;
    adupsoil_layer_day_point = adup_noPALMS;
    load('noPALMSsoy');
    grain_carbon_day_point = grainC_soy./1000;
    leaf_carbon_day_point = leafC_soy./1000;
    rhizome_carbon_day_point = rhizC_soy./1000;
    root_carbon_day_point = rootC_soy./1000;
    stem_carbon_day_point = stemC_soy./1000;
end


%     
% % From PALMS, the units are / smp [MPa] /layeruptake [mm/d] / qq [m/d]
% for i=start_year:end_year
%     a='data2NcycleLighten';
%     b='dataCarbon2NcycleLighten';
%     c=num2str(i);
%     aa='_S';
%     bb=num2str(scenarios);
%     d=strcat(a,c,aa,bb);
%     e=strcat(b,c,aa,bb);
%     
%     % load PALMS' output except for Carbon contents
%     for k=4 %1:size(Need2becalled,2) %susana edit SOIL MOISTURE CALLED
%         f=Need2becalled(k);
%         g=num2str(choose_point);
%         h=strcat(f,g);
%         p=char(h);
%         o=load(d,p);
%         
%         q=eval(['o' '.' p]);
%         % totlail is not a metrix, but a vector
%         if k==1
%             q=q(1,:);
%         end
%         % delete 26th layers due to the number of water flux between
%         % layers(26 layers)
%         if k>2
%             q(26,:)=[];
%         end
%         % Put the year at the end of loaded variables
%         eval([p '_' num2str(i) '= q;']);
%         
%     end
% %     % load only PALMS' Carbon contents
% %     for k=1:size(Need2becalledCarbon,2)
% %         f=Need2becalledCarbon(k);
% %         g=num2str(choose_point);
% %         h=strcat(f,g);
% %         p=char(h);
% %         o=load(e,p);
% %         
% %         % Put the year at the end of loaded variables
% %         q=eval(['o' '.' p]);
% %         % Carbon contents have the specific plant funcation types
% %         if choose_crop==2
% %             for j=1:size(soy_year,2)
% %                 if i == soy_year(j)
% %                     choose_pft=choose_pft_Soy;
% %                 end
% %             end
% %             for j=1:size(corn_year,2)
% %                 if i == corn_year(j)
% %                     choose_pft=choose_pft_Corn;
% %                 end
% %             end
% %         end
% %         q=q(choose_pft,:);
% %         eval([p '_' num2str(i) '= q;']);
% %     end
% end
% % disp('Done for loading all the needed data!')
% 
% %%-----------------------------------------------------------------------%%
% % 3. Putting needed variables into single matrixs
% % for 'Need2becalled' + 'Need2becalledCarbon' + TBMLa
% %CombinedNeed2becalled=[Need2becalled, Need2becalledCarbon, 'TBMLa'];
% CombinedNeed2becalled=[Need2becalled, Need2becalledCarbon];
% 
% disp('Organizing the data')
% for k=4 %1:size(CombinedNeed2becalled,2) SOIL MOISTURE ONLY
    % Pre-allocate Matrix that will be used
%     for i=start_year:end_year
%         f=CombinedNeed2becalled(k);
%         g=num2str(choose_point);
%         h='_';
%         o=num2str(i);
%         p=strcat(f,g,h,o);
%         q=char(p);
%         s=char(f);
%         r=eval([q]);
%         
%         Size_This=size(r,1);
%         eval( [s '= zeros(Size_This);'] );
%     end
%     
    % Check point
%     count_year=0;
%     for i=start_year:end_year
%         count_year=count_year+1;
%         
%         f=CombinedNeed2becalled(k);
%         g=num2str(choose_point);
%         h='_';
%         o=num2str(i);
%         p=strcat(f,g,h,o);
%         q=char(p);
%         
%         r=eval([q]);
%         % Check point
%         s=char(f);
%         
%         for j=1:365
%             eval( [s '(:,(count_year-1)*365+j)'  '= r(:,j);'] );
%         end
%     end
% end
% disp('Done for organizing data')
%END SUSANA EDIT

%%-----------------------------------------------------------------------%%
% 4. Load weather forcings. Especially precipitation & wind speed, which
% will be used in volatilization
load('Daily_PPT_1798_2001')
load('Daily_Wind_1798_2001')
load('Daily_PPT_2002_2011')
load('Daily_Wind_2002_2011')
if scenarios == 1
    load('Daily_PPT_2012_2012_S1')
    load('Daily_Wind_2012_2012_S1')
    
    Ppt_1798_2012_mm=[Daily_PPT_1798_2001 Daily_PPT_2002_2011 Daily_PPT_2012_2012];
    Wind_1798_2012_ms=[Daily_Wind_1798_2001 Daily_Wind_2002_2011 Daily_Wind_2012_2012];
    
    PPT_mm=Ppt_1798_2012_mm((start_year-1798)*365+1:size(Ppt_1798_2012_mm,2)-(2012-end_year)*365);
    Wind_ms=Wind_1798_2012_ms((start_year-1798)*365+1:size(Wind_1798_2012_ms,2)-(2012-end_year)*365);
elseif scenarios == 2
    load('Daily_PPT_2012_2012_S2')
    load('Daily_Wind_2012_2012_S2')
    
    Ppt_1798_2012_mm=[Daily_PPT_1798_2001 Daily_PPT_2002_2011 Daily_PPT_2012_2012];
    Wind_1798_2012_ms=[Daily_Wind_1798_2001 Daily_Wind_2002_2011 Daily_Wind_2012_2012];
    
    PPT_mm=Ppt_1798_2012_mm((start_year-1798)*365+1:size(Ppt_1798_2012_mm,2)-(2012-end_year)*365);
    Wind_ms=Wind_1798_2012_ms((start_year-1798)*365+1:size(Wind_1798_2012_ms,2)-(2012-end_year)*365);
elseif scenarios == 3
    load('Daily_PPT_2012_2012_S3')
    load('Daily_Wind_2012_2012_S3')
    Ppt_1798_2012_mm=[Daily_PPT_1798_2001 Daily_PPT_2002_2011 Daily_PPT_2012_2012];
    Wind_1798_2012_ms=[Daily_Wind_1798_2001 Daily_Wind_2002_2011 Daily_Wind_2012_2012];
    
    PPT_mm=Ppt_1798_2012_mm((start_year-1798)*365+1:size(Ppt_1798_2012_mm,2)-(2012-end_year)*365);
    Wind_ms=Wind_1798_2012_ms((start_year-1798)*365+1:size(Wind_1798_2012_ms,2)-(2012-end_year)*365);
elseif scenarios == 4
    load('Daily_PPT_2012_2012_S4')
    load('Daily_Wind_2012_2012_S4')
    Ppt_1798_2012_mm=[Daily_PPT_1798_2001 Daily_PPT_2002_2011 Daily_PPT_2012_2012];
    Wind_1798_2012_ms=[Daily_Wind_1798_2001 Daily_Wind_2002_2011 Daily_Wind_2012_2012];
    
    PPT_mm=Ppt_1798_2012_mm((start_year-1798)*365+1:size(Ppt_1798_2012_mm,2)-(2012-end_year)*365);
    Wind_ms=Wind_1798_2012_ms((start_year-1798)*365+1:size(Wind_1798_2012_ms,2)-(2012-end_year)*365);
elseif scenarios == 5
    load('Daily_PPT_2012_2050_S5')
    load('Daily_Wind_2012_2050_S5')
    load('Daily_PPT_2051_2107_S5')
    load('Daily_Wind_2051_2107_S5')
    Ppt_1798_2107_mm=[Daily_PPT_1798_2001 Daily_PPT_2002_2011 Daily_PPT_2012_2050 Daily_PPT_2051_2107];
    Wind_1798_2107_ms=[Daily_Wind_1798_2001 Daily_Wind_2002_2011 Daily_Wind_2012_2050 Daily_Wind_2051_2107];
    
    PPT_mm=Ppt_1798_2107_mm((start_year-1798)*365+1:size(Ppt_1798_2107_mm,2)-(2107-end_year)*365);
    Wind_ms=Wind_1798_2107_ms((start_year-1798)*365+1:size(Wind_1798_2107_ms,2)-(2107-end_year)*365);    
end



%%-----------------------------------------------------------------------%%
% 5. Calculate biomass drops by Woo et al., 2014.
CN_TotalBomassDrop_CNratio();
    if Fix_Crop_Vaues==1
        TBMLa=pre_TBMLa;
        TBMLaNoLitterBack=pre_TBMLaNoLitterBack;
        RootDeathWhenHarvest=pre_RootDeathWhenHarvest;
        TBMLaOnlyLitterBack=pre_TBMLaOnlyLitterBack; 
        TBMLaNoLitterBack=pre_TBMLaNoLitterBack; 
    end



FORCING.TBMla=TBMLa;
FORCING.TBMLaOnlyLitterBack=TBMLaOnlyLitterBack; % Only harvest litter back to soil
FORCING.TBMLaNoLitterBack=TBMLaNoLitterBack; % Only natural litter back to soil

%%-----------------------------------------------------------------------%%
% 6. Put C/N ratio in biomass pool
for i=1:PARAMS.nl_soil
    CNb(i)=CN_biomass;
end
CNb=CNb(:);
PARAMS.CN.CNb=CNb;

%%-----------------------------------------------------------------------%%
% 7. Compute total carobn in above- and belowground
Total_Above_Carbon=grain_carbon_day_point+leaf_carbon_day_point+...
    stem_carbon_day_point;
Total_Below_Carbon=rhizome_carbon_day_point+root_carbon_day_point;

%%-----------------------------------------------------------------------%%
% 8. Initial value. It does not matter what you put, but matter of time to
% reach the steady state. 
if from_previous_data ==0
    % Now, we can put VARIABLES.Nl, VARIABLES.Nh and VARIABLES.CNl.
    %SUSANA EDIT - ADJUSTED FOR 12 SOIL LAYERS (REMOVE PALMS)
%     VARIABLES.Nl=[195.4452;140.0438;117.1134;107.8897;106.2631;...
%         106.3852;106.6094;102.0573;96.7170;92.5914;86.8880;...
%         73.3561;57.6097;42.8340;33.6560;30.7472;24.9178;16.3360;...
%         11.0637;7.5685;5.1329;3.2966;2.0030;0.5694;0.0360];
    VARIABLES.N1=[167.7445;112.50155;106.3;104.33335;94.6542;78.9;50.22185;32.2016;...
        20.6269;9.3161;4.1;0.85612];

    %VARIABLES.Nl=VARIABLES.Cl./20;
    %SUSANA EDIT - ADJUSTED FOR 12 SOIL LAYERS (REMOVE PALMS)
%     VARIABLES.Nh=[926.7950;925.4775;924.1605;922.8435;921.6120;...
%         905.1880;875.4310;835.4210;787.0845;732.3465;660.9225;...
%         534.4745;406.6825;297.3605;221.9140;195.7475;154.7765;...
%         103.0470;68.6060;45.6765;30.4105;19.2425;11.5725;3.2460;0.4245];
    VARIABLES.Nh=[926.13625;923.502;912.5;855.426;759.7155;586.5;352.0215;208.83075;...
        128.91175;57.14125;24.2;4.9113];

    if SWITCHES.CN_type == 0
        VARIABLES.Nh=nan(nl_soil,1);
    end
    %SUSANA EDIT - ONLY 12 LAYERS (REMOVE PALMS)
%     VARIABLES.CNl=[9.8112;21.2002;18.1108;17.7465;17.7009;17.7192;...
%         17.8363;17.8280;17.8546;17.9807;18.1437;18.3011;18.4070;...
%         18.4442;18.6144;18.7797;18.9235;18.9192;19.0357;19.1923;...
%         19.3245;19.4537;19.5885;19.6925;13.4932];
    VARIABLES.CN1=[15.5057;17.92865;17.7;17.83215;17.91765;18.2;18.4256;...
        18.69705;18.92135;19.114;19.4;19.6717];


    %VARIABLES.CNl=ones(size(VARIABLES.Nl,1),1)*20;
end

%%-----------------------------------------------------------------------%%
% 9. Now is good time to compute kl, kd, and kh if you do not have the
% previous simulation values. FYI, kl, kd, and kh are decomposition
% parameters and you can find the euqations in Porporato et al., 2003. 
if from_previous_data ==0
    % Calculate kl, kd, and kh with corn and soy
    [kl, kd, kh] = CN_Compute_kl_kd_kh (PARAMS, VARIABLES, FORCING...
        , VERTSTRUC, SWITCHES, CONSTANTS, vwc_layer_day_point, ice_layer_day_point, temp_layer_day_point...
        , TBMLa, rootfr, WeithedAveCNabove, WeithedAveCNbelow, WeithedAveCNcrRatio, RootDeathWhenH, soiltype, choose_crop);
    
    % Negative? kl, kh, kd?
    disp('please check Negative kl, kh, and kd')
    
    PARAMS.CN.kl=kl;
    PARAMS.CN.kd=kd;
    PARAMS.CN.kh=kh;
    
    % Set reasonable Amm and NH3 amount
    [Amm, NH3, Amm_mo, Amm_im, VARIABLES] = CN_compute_initial_Amm_NH3(PARAMS, temp_layer_day_point, nl_soil, VARIABLES);
    
   
end

%%-----------------------------------------------------------------------%%
% 9. Multiply necessary varibles for initialting Cl, Cb Ch, Nit, and Amm
% I.e., to get the steady state. Please check this. You can only use this
% when you do not have enough forcing to simulate.
How_many_times_to_repmat=ceil(how_many_year/(end_year-start_year+1));

totlail_layer_day_point=repmat(totlail_layer_day_point,1,How_many_times_to_repmat);
LAI_All_layer=repmat(totlail_layer_day_point,nl_soil,How_many_times_to_repmat);

vwc_layer_day_point=repmat(vwc_layer_day_point,1,How_many_times_to_repmat);
ice_layer_day_point=repmat(ice_layer_day_point,1,How_many_times_to_repmat);
temp_layer_day_point=repmat(temp_layer_day_point,1,How_many_times_to_repmat);
adupsoil_layer_day_point=repmat(adupsoil_layer_day_point,1,How_many_times_to_repmat);
qq_layer_day_point=repmat(qq_layer_day_point,1,How_many_times_to_repmat);

grain_carbon_day_point=repmat(grain_carbon_day_point,1,How_many_times_to_repmat);
leaf_carbon_day_point=repmat(leaf_carbon_day_point,1,How_many_times_to_repmat);
rhizome_carbon_day_point=repmat(rhizome_carbon_day_point,1,How_many_times_to_repmat);
root_carbon_day_point=repmat(root_carbon_day_point,1,How_many_times_to_repmat);
stem_carbon_day_point=repmat(stem_carbon_day_point,1,How_many_times_to_repmat);

TBMLa=repmat(TBMLa,1,How_many_times_to_repmat);
CNveg=repmat(CNratio_apply,1,How_many_times_to_repmat);
CNveg_root=repmat(CNratio_root,1,How_many_times_to_repmat);
if choose_crop==3
    load('MGSG_BelowCN','SG_Below_CN')
    CNratio_root=repmat((SG_Below_CN-10)',1,how_many_year*How_many_times_to_repmat);
    CNveg_root=repmat((SG_Below_CN-10)',1,how_many_year*How_many_times_to_repmat);
end
if choose_crop==4
    load('MGSG_BelowCN','MG_Below_CN')
    CNratio_root=repmat((MG_Below_CN)',1,how_many_year*How_many_times_to_repmat);
    CNveg_root=repmat((MG_Below_CN)',1,how_many_year*How_many_times_to_repmat);
end

Drain_NP=repmat(Drain_NP,1,How_many_times_to_repmat);
Drain_Corn=repmat(Drain_Corn,1,How_many_times_to_repmat);
Drain_MG=repmat(Drain_MG,1,How_many_times_to_repmat);
Drain_SG=repmat(Drain_SG,1,How_many_times_to_repmat);

if choose_crop==1
    Drain=Drain_NP;
elseif choose_crop==2
    Drain=Drain_Corn;
elseif choose_crop==3
    Drain=Drain_SG;
elseif choose_crop==4
    Drain=Drain_MG;
end
    
    
%%-----------------------------------------------------------------------%%
% 10. Put the forcings in structure
FORCING.TBMla=TBMLa;

% CNveg will be used below, so change the name
ALL_CNveg=CNveg;
clear CNveg
ALL_CNveg_root=CNveg_root;
clear CNveg_root

% save carbon variables in structure
if Fix_Crop_Vaues==1
    grain_carbon_day_point=pre_grain_carbon_day_point;
    leaf_carbon_day_point=pre_leaf_carbon_day_point;
    rhizome_carbon_day_point=pre_rhizome_carbon_day_point;
    root_carbon_day_point=pre_root_carbon_day_point;
    stem_carbon_day_point=pre_stem_carbon_day_point;
end

CARBONS.grain_carbon_day_point=grain_carbon_day_point;
CARBONS.leaf_carbon_day_point=leaf_carbon_day_point;
CARBONS.rhizome_carbon_day_point=rhizome_carbon_day_point;
CARBONS.root_carbon_day_point=root_carbon_day_point;
CARBONS.stem_carbon_day_point=stem_carbon_day_point;
clear grain_carbon_day_point
clear leaf_carbon_day_point
clear rhizome_carbon_day_point
clear root_carbon_day_point
clear stem_carbon_day_point
Total_Carbon=CARBONS.grain_carbon_day_point+CARBONS.leaf_carbon_day_point+...
    CARBONS.rhizome_carbon_day_point+CARBONS.root_carbon_day_point+...
    CARBONS.stem_carbon_day_point;

%%-----------------------------------------------------------------------%%
% 11. Fertilizer
% Preallocate the fertilzer inputs
Nit_fertilizer1=zeros(1,365);
Nit_fertilizer2=zeros(1,365);
Nit_fertilizer3=zeros(1,365);
Nit_fertilizer4=zeros(1,365);
Nit_fertilizer5=zeros(1,365);
Nit_fertilizer6=zeros(1,365);

Amm_fertilizer1=zeros(1,365);
Amm_fertilizer2=zeros(1,365);
Amm_fertilizer3=zeros(1,365);
Amm_fertilizer4=zeros(1,365);
Amm_fertilizer5=zeros(1,365);
Amm_fertilizer6=zeros(1,365);

Organic_N_fertilizer1=zeros(1,365);
Organic_N_fertilizer2=zeros(1,365);
Organic_N_fertilizer3=zeros(1,365);
Organic_N_fertilizer4=zeros(1,365);
Organic_N_fertilizer5=zeros(1,365);
Organic_N_fertilizer6=zeros(1,365);

% Note that the unit is g/m2
Nit_fertilizer1(Day_of_fertilzer_1st_year)=Nitrate_fertilizer_amount_1st_year;
Nit_fertilizer2(Day_of_fertilzer_2nd_year)=Nitrate_fertilizer_amount_2nd_year;
Nit_fertilizer3(Day_of_fertilzer_3rd_year)=Nitrate_fertilizer_amount_3rd_year;
Nit_fertilizer4(Day_of_fertilzer_4th_year)=Nitrate_fertilizer_amount_4th_year;
Nit_fertilizer5(Day_of_fertilzer_5th_year)=Nitrate_fertilizer_amount_5th_year;
Nit_fertilizer6(Day_of_fertilzer_6th_year)=Nitrate_fertilizer_amount_6th_year;

Amm_fertilizer1(Day_of_fertilzer_1st_year)=Ammonium_fertilizer_amount_1st_year;
Amm_fertilizer2(Day_of_fertilzer_2nd_year)=Ammonium_fertilizer_amount_2nd_year;
Amm_fertilizer3(Day_of_fertilzer_3rd_year)=Ammonium_fertilizer_amount_3rd_year;
Amm_fertilizer4(Day_of_fertilzer_4th_year)=Ammonium_fertilizer_amount_4th_year;
Amm_fertilizer5(Day_of_fertilzer_5th_year)=Ammonium_fertilizer_amount_5th_year;
Amm_fertilizer6(Day_of_fertilzer_6th_year)=Ammonium_fertilizer_amount_6th_year;

Organic_N_fertilizer1(Day_of_fertilzer_1st_year)=Organic_N_fertilizer_1st_year;
Organic_N_fertilizer2(Day_of_fertilzer_2nd_year)=Organic_N_fertilizer_2nd_year;
Organic_N_fertilizer3(Day_of_fertilzer_3rd_year)=Organic_N_fertilizer_3rd_year;
Organic_N_fertilizer4(Day_of_fertilzer_4th_year)=Organic_N_fertilizer_4th_year;
Organic_N_fertilizer5(Day_of_fertilzer_5th_year)=Organic_N_fertilizer_5th_year;
Organic_N_fertilizer6(Day_of_fertilzer_6th_year)=Organic_N_fertilizer_6th_year;

% 1,2,3 for current use, 4,5,6 for future
Input_of_Nit_fertilizer1=[Nit_fertilizer1 Nit_fertilizer2 Nit_fertilizer3];
Input_of_Amm_fertilizer1=[Amm_fertilizer1 Amm_fertilizer2 Amm_fertilizer3];
Input_of_Organic_fertilizer1=[Organic_N_fertilizer1 Organic_N_fertilizer2 Organic_N_fertilizer3];
Input_of_Nit_fertilizer2=[Nit_fertilizer4 Nit_fertilizer5 Nit_fertilizer6];
Input_of_Amm_fertilizer2=[Amm_fertilizer4 Amm_fertilizer5 Amm_fertilizer6];
Input_of_Organic_fertilizer2=[Organic_N_fertilizer4 Organic_N_fertilizer5 Organic_N_fertilizer6];

% put Fertilizer as much as the number of days
Repmat_N=floor(how_many_year/3);
IN_Nit_Fertilizer=repmat(Input_of_Nit_fertilizer2,1,Repmat_N);
IN_Amm_Fertilizer=repmat(Input_of_Amm_fertilizer2,1,Repmat_N);
IN_Organic_Fertilizer=repmat(Input_of_Organic_fertilizer2,1,Repmat_N);

The_rest=(how_many_year/3)-Repmat_N;
if The_rest>0.333
    IN_Nit_Fertilizer=[IN_Nit_Fertilizer Nit_fertilizer4];
    IN_Amm_Fertilizer=[IN_Amm_Fertilizer Amm_fertilizer4];
    IN_Organic_Fertilizer=[IN_Organic_Fertilizer Organic_N_fertilizer4];
    
    if The_rest>0.666
        IN_Nit_Fertilizer=[IN_Nit_Fertilizer Nit_fertilizer5];
        IN_Amm_Fertilizer=[IN_Amm_Fertilizer Amm_fertilizer5];
        IN_Organic_Fertilizer=[IN_Organic_Fertilizer Organic_N_fertilizer5];
    end
end
IN_Nit_Fertilizer=[Input_of_Nit_fertilizer1 IN_Nit_Fertilizer];
IN_Amm_Fertilizer=[Input_of_Amm_fertilizer1 IN_Amm_Fertilizer];
IN_Organic_Fertilizer=[Input_of_Organic_fertilizer1 IN_Organic_Fertilizer];

% This is only for corn-soybean roation. For EBI site, corn-soybean had
% been cultivated.
if initialization == 1
    Repmat_NN=ceil(how_many_year/2);
    
    Nit_fertilizer_Initial=[Nit_fertilizer5 Nit_fertilizer6];
    Nit_fertilizer_Initial=repmat(Nit_fertilizer_Initial,1,Repmat_NN);
    IN_Nit_Fertilizer=[Nit_fertilizer_Initial Nit_fertilizer1 Nit_fertilizer2 Nit_fertilizer3 Nit_fertilizer4];
    
    Amm_fertilizer_Initial=[Amm_fertilizer5 Amm_fertilizer6];
    Amm_fertilizer_Initial=repmat(Amm_fertilizer_Initial,1,Repmat_NN);
    IN_Amm_Fertilizer=[Amm_fertilizer_Initial Amm_fertilizer1 Amm_fertilizer2 Amm_fertilizer3 Amm_fertilizer4];
    
    Organic_fertilizer_Initial=[Organic_N_fertilizer5 Organic_N_fertilizer6];
    Organic_fertilizer_Initial=repmat(Organic_fertilizer_Initial,1,Repmat_NN);
    IN_Organic_Fertilizer=[Organic_fertilizer_Initial Organic_N_fertilizer1 Organic_N_fertilizer2 Organic_N_fertilizer3 Organic_N_fertilizer4];
end

%%-----------------------------------------------------------------------%%
% 12. Fixation 
% repmat for Soy and switchgrass (zeros)
nfix1=ones(1,365)*nfix(1);
nfix2=ones(1,365)*nfix(2);
nfix3=ones(1,365)*nfix(3);

In_nfix1=[nfix1 nfix2 nfix3];
In_nfix2=repmat(In_nfix1,1,Repmat_N);
if The_rest>0.333
    nfix_rate=[In_nfix2 nfix1];
    if The_rest>0.666
        nfix_rate=[nfix_rate nfix2];
    end
else
    nfix_rate=In_nfix2;
end
PARAMS.nfix_rate=nfix_rate;

% for miscanthus
if choose_crop == 4
    nfix4=ones(1,365)*nfix(4);
    nfix5=ones(1,365*(how_many_year-4))*nfix(5);
    PARAMS.nfix_rate=[nfix1 nfix2 nfix3 nfix4 nfix5];
end

% for initialization, corn-soybean rotation
if initialization == 1
    nfix_initial=[nfix2 nfix3];
    nfix_initial=repmat(nfix_initial,1,Repmat_NN);
    nfix_rate=[nfix_initial nfix2 nfix2 nfix3 nfix2];
    
    PARAMS.nfix_rate=nfix_rate;
end

%%-----------------------------------------------------------------------%%
% 13. N remobilization
if SWITCHES.Nremobil
    Nremobilization1=ones(1,365)*Nremobilization(1);
    Nremobilization2=ones(1,365)*Nremobilization(2);
    Nremobilization3=ones(1,365)*Nremobilization(3);
    Nremobilization4=ones(1,365)*Nremobilization(4);
    Nremobilization5=ones(1,365*(how_many_year-4))*Nremobilization(5);
    PARAMS.Nremobilization_rate=[Nremobilization1 Nremobilization2 Nremobilization3 Nremobilization4 Nremobilization5];
end

%%-----------------------------------------------------------------------%%
% 14. Put atmospheric N deposition at many as the number of days
% 1,2,3 4 for current use(08 09 10 11), 5 for future [g/m2/d]
if SWITCHES.CN.Ndeposition == 1
    for i=1:365
        YDeposition_NH4_N1(i)=Deposition_NH4_N(1);
        YDeposition_NH4_N2(i)=Deposition_NH4_N(2);
        YDeposition_NH4_N3(i)=Deposition_NH4_N(3);
        YDeposition_NH4_N4(i)=Deposition_NH4_N(4);
        YDeposition_NH4_N5(i)=Deposition_NH4_N(5);
        
        YDeposition_NO3_N1(i)=Deposition_NO3_N(1);
        YDeposition_NO3_N2(i)=Deposition_NO3_N(2);
        YDeposition_NO3_N3(i)=Deposition_NO3_N(3);
        YDeposition_NO3_N4(i)=Deposition_NO3_N(4);
        YDeposition_NO3_N5(i)=Deposition_NO3_N(5);
    end
    Deposition_NH4_N_Apply1=repmat(YDeposition_NH4_N5,1,how_many_year);
    Deposition_NO3_N_Apply1=repmat(YDeposition_NO3_N5,1,how_many_year);
end
Deposition_NH4_N_Apply=[YDeposition_NH4_N1 YDeposition_NH4_N2...
    YDeposition_NH4_N3 YDeposition_NH4_N4 Deposition_NH4_N_Apply1];
Deposition_NO3_N_Apply=[YDeposition_NO3_N1 YDeposition_NO3_N2...
    YDeposition_NO3_N3 YDeposition_NO3_N4 Deposition_NO3_N_Apply1];

if initialization == 1
    Deposition_NH4_N_Apply=[Deposition_NH4_N_Apply1 YDeposition_NH4_N1 YDeposition_NH4_N2...
        YDeposition_NH4_N3 YDeposition_NH4_N4];
    Deposition_NO3_N_Apply=[Deposition_NO3_N_Apply1 YDeposition_NO3_N1 YDeposition_NO3_N2...
        YDeposition_NO3_N3 YDeposition_NO3_N4];
end


%%-----------------------------------------------------------------------%%
% 16. Put the previous data. If you have the previous data. 
if from_previous_data == 1
    VARIABLES.Cl=Previous_Cl;
    Cl=Previous_Cl;
    VARIABLES.Cb=Previous_Cb;
    Cb=Previous_Cb;
    VARIABLES.Ch=Previous_Ch;
    Ch=Previous_Ch;
    VARIABLES.CNl=Previous_CNl;
    CNl=Previous_CNl;
    VARIABLES.Nl=Previous_Nl;
    Nl=Previous_Nl;
    
    %------------------------------
    % Need chaging
    VARIABLES.Amm=Previous_Amm;
    Amm=Previous_Amm;

    
    VARIABLES.Amm_mo=Previous_Amm_mo;
    Amm_mo=VARIABLES.Amm_mo;
        
    VARIABLES.Amm_im=Previous_Amm_im;
    Amm_im=VARIABLES.Amm_im;
        
    Amm_mo2im = zeros(nl_soil,1);
    Amm_im2mo = zeros(nl_soil,1);
    %------------------------------
    
    VARIABLES.Nit=Previous_Nit;
    Nit=Previous_Nit;
    
    VARIABLES.NH3=Previous_NH3;
    NH3=Previous_NH3;
    %VARIABLES.pH=6.*ones(nl_soil,1);
    %pH=VARIABLES;
    % Until here!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    PARAMS.CN.kl=Previous_kl;
    kl=Previous_kl;
    PARAMS.CN.kd=Previous_kd;
    kd=Previous_kd;
    PARAMS.CN.kh=Previous_kh;
    kh=Previous_kh;
    
    VARIABLES.Nh=VARIABLES.Ch./PARAMS.CN.CNh;
end

%%-----------------------------------------------------------------------%%
% 17. If you would like to fix kl, kd, and kh, decomposition parameters, 
% from the previous simulation.
if fixed_kl_kd_kh == 1;
    PARAMS.CN.kl=Fixed_kl;
    kl=Fixed_kl;
    PARAMS.CN.kd=Fixed_kd;
    kd=Fixed_kd;
    PARAMS.CN.kh=Fixed_kh;
    kh=Fixed_kh;
end