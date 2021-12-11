%*************************************************************************
%               Compute Biomass Drop in Carbon & CN ratio
%*************************************************************************

%%-----------------------------------------------------------------------%%
% 1. Compute the first day of year
for i=1:size(leaf_carbon_day_point,2)-1
    % Compute Year
    for ii=1:1:how_many_year
        BD_detect_first_day_of_year(1)=1;
        BD_detect_first_day_of_year(ii+1)=1+365*ii;
    end
    BD_detect_first_day_of_year(end)=[];
end

%%-----------------------------------------------------------------------%%
% 2. Calculate harvest and emergence day
for i=1:size(leaf_carbon_day_point,2)/365
    k1=365*(i-1)+1;
    
    % Check Harvest day
    current_year_LC=leaf_carbon_day_point(k1:k1+364);
    FirstTime=0;
    for j=200:365
        if FirstTime == 0
            if current_year_LC(j) < 0.001
                FirstTime=FirstTime+1;
                harvest_day=j;
            end
        end
    end
    
    if i == 1
        HD_within_year=harvest_day;
        HD_over_year=harvest_day;
    else
        HD_within_year=[HD_within_year harvest_day];
        HD_over_year=[HD_over_year harvest_day+365*(i-1)];
    end
    
    % Check Emergence day
    FirstTimeE=0;
    for j=1:365
        if FirstTimeE == 0
            if current_year_LC(j) > current_year_LC(1)
                FirstTimeE=FirstTimeE+1;
                Emergence_day=j;
            end
        end
    end
    
    if i == 1
        ED_within_year=Emergence_day;
        ED_over_year=Emergence_day;
    else
        ED_within_year=[ED_within_year Emergence_day];
        ED_over_year=[ED_over_year Emergence_day+365*(i-1)];
    end
end

%%-----------------------------------------------------------------------%%
% 3. For corn and soybean grain carbon is not considered as source for the
% harvest residue drop due to harvest practice.
dtime=86400;
BD_Check_year=0;

if choose_crop==2
    Tot_AboveC=(leaf_carbon_day_point+stem_carbon_day_point)*1000; % [kg/m2]
end
if choose_crop==3
    Tot_AboveC=(leaf_carbon_day_point+stem_carbon_day_point...
        +grain_carbon_day_point)*1000; % [kg/m2]
end
if choose_crop==4
    Tot_AboveC=(leaf_carbon_day_point+stem_carbon_day_point...
        +grain_carbon_day_point)*1000; % [kg/m2]
end
Tot_BelowC=(root_carbon_day_point+rhizome_carbon_day_point)*1000; % [kg/m2]

%%-----------------------------------------------------------------------%%
% 4. Applying shift biomass drop!: litterfall (Woo et al., 2014) 
BD_Year=start_year-1;
for i=1:size(leaf_carbon_day_point,2)/365
    BD_Year=BD_Year+1;
    
    if choose_crop==2
        for j=1:size(soy_year,2)
            if BD_Year == soy_year(j)
                leafspan=leafspan_Soy;
                HavestResidueRate=HavestResidueRate_Soy;
            end
        end
        for j=1:size(corn_year,2)
            if BD_Year == corn_year(j)
                leafspan=leafspan_Corn;
                HavestResidueRate=HavestResidueRate_Corn;
            end
        end
    end

    k1=365*(i-1)+1;
    
    current_year_LC=leaf_carbon_day_point(k1:k1+364).*1000;
    ShiftLC(1+leafspan:365+leafspan)=current_year_LC;
    
    counthere=0;
    for j=ED_within_year(i):365
        Pre_BD=ShiftLC(j+1)-ShiftLC(j);
        if j == 280
            pause = 1;
        end
        if Pre_BD > 0
            counthere=counthere+1;
            BD(j)=Pre_BD;
            if counthere == 1
                BD(j)=0;
            end
        else
            BD(j)=0;
        end
    end
    
    BD(HD_within_year(i):365)=0;
    
    if i == 1
        Bdrop=BD;
    else
        Bdrop=[Bdrop BD];
    end
    
end
TBMLaNoLitterBack=Bdrop;

%%-----------------------------------------------------------------------%%
% 5. Harvest biomass drop (Woo et al., 2014) 
for i=1:size(leaf_carbon_day_point,2)-1
    for j = 1:size(HD_over_year,2)
        if i  == HD_over_year(j)
            % TBMLa only Litter back when harvest
            TBMLaOnlyLitterBack(i)=Tot_AboveC(i).*HavestResidueRate;
            TBMLaFromLeaf(i) = leaf_carbon_day_point(i).*1000.*HavestResidueRate;
            TBMLaFromStem(i) = stem_carbon_day_point(i).*1000.*HavestResidueRate;
            TBMLaFromGrain(i) = 0;
            
            if choose_crop ~= 2
                TBMLaFromGrain(i) = grain_carbon_day_point(i).*1000.*HavestResidueRate;
            end
             
        end
    end
end

% For the gap
GapDay=size(Bdrop,2)-size(TBMLaOnlyLitterBack,2);
TBMLaOnlyLitterBack=[TBMLaOnlyLitterBack zeros(1,GapDay)];
TBMLaFromLeaf=[TBMLaFromLeaf zeros(1,GapDay)];
TBMLaFromStem=[TBMLaFromStem zeros(1,GapDay)];
TBMLaFromGrain=[TBMLaFromGrain zeros(1,GapDay)];

%%-----------------------------------------------------------------------%%
% 6. Total aboveground litter input into the soil 
TBMLa=Bdrop+TBMLaOnlyLitterBack;

% Check the cumulative drops.
for i=1:size(BD_detect_first_day_of_year,2)
    CumulBdropCheck=cumsum(TBMLa(365*(i-1)+1:365*i));
    CumulBdropEachYearCheck(365*(i-1)+1:365*i)=CumulBdropCheck;
end

%%-----------------------------------------------------------------------%%
% 7. Harvest belowgroudn residue drop due to harvest practice.
RootDeathWhenHarvest=zeros(1,size(leaf_carbon_day_point,2));
CountRoot=0;
for i=1:size(leaf_carbon_day_point,2)
    for j = 1:size(HD_over_year,2)
        if i  == HD_over_year(j)
            % Put the litter back when it is harvested
            %TBMLa(i)=Tot_AboveC(i).*HavestResidueRate;  
            %TBMLaOnlyLitterBack(i)=Tot_AboveC(i).*HavestResidueRate;            
            
            % Compute root death for Corn and Soy
            if choose_crop==2
                CountRoot=CountRoot+1;
                
                RootDeathWhenHarvest(i)=Tot_BelowC(i);                
                RootDeathWhenH(CountRoot)=Tot_BelowC(i);    
            end
            if choose_crop==3
                CountRoot=CountRoot+1;
                
                RootDeathWhenHarvest(i)=0;                
                RootDeathWhenH(CountRoot)=0;    
            end
            if choose_crop==4
                CountRoot=CountRoot+1;
                
                RootDeathWhenHarvest(i)=0;                
                RootDeathWhenH(CountRoot)=0;    
            end
        end
    end
end

%%-----------------------------------------------------------------------%%
% 8. Compute C/N ratio
% For applying C:N ratio
BD_j_day=0;
BD_Check_year=0;
for i=1:size(leaf_carbon_day_point,2)
    % Compute Year
    for ii=1:1:how_many_year
        BD_detect_first_day_of_year(1)=1;
        BD_detect_first_day_of_year(ii+1)=1+365*ii;
    end
    BD_detect_first_day_of_year(end)=[];
    for ii=1:how_many_year
        if i == BD_detect_first_day_of_year(ii)
            BD_Check_year=BD_Check_year+1;
            % Check How many years has been run!
            BD_Year=start_year+BD_Check_year-1;
            BD_j_day=0;
        end
    end
    % Each Year Day
    BD_j_day=BD_j_day+1;
    
    % Corn-corn-soybean
    if choose_crop==2
        for j=1:size(soy_year,2)
            if BD_Year == soy_year(j)
                CN_leaf=CN_leaf_Soy;
                CN_stem=CN_stem_Soy;
                CN_grain=CN_grain_Soy;
                CN_root=CN_root_Soy;
                
                CR_ratio1=PARAMS.CR_ratio_Soy;
            end
        end
        for j=1:size(corn_year,2)
            if BD_Year == corn_year(j)
                CN_leaf=CN_leaf_Corn;
                CN_stem=CN_stem_Corn;
                CN_grain=CN_grain_Corn;
                CN_root=CN_root_Corn;
                
                CR_ratio1=PARAMS.CR_ratio_Corn;
            end
        end
    else
        CR_ratio1=PARAMS.CN.CR_ratio;
    end
    
    % Switchgrass
    if choose_crop==3
        CN_leaf=21.154*exp(0.0023*BD_j_day); % Original - low MIN
        
        %CN_stem=0.691*BD_j_day-65.603; % Linear
        CN_stem=0.8375*BD_j_day-121.95; % Linear
        if CN_stem < 5
            CN_stem = 5;
        end
        CN_grain=CN_stem;
    end
    
    % Miscanthus
    if choose_crop==4
        %CN_leaf=13.525*exp(0.0046*BD_j_day); % Original - low MIN
        CN_leaf=15.522*exp(0.0035*BD_j_day); % GOOD
        %CN_leaf=14.333*exp(0.0039*BD_j_day);
        %CN_leaf=16.031*exp(0.0033*BD_j_day);
        
        %CN_stem=0.4059*exp(0.0193*BD_j_day); % Exponential
        CN_stem=1.088*BD_j_day-171.91; % Linear
        if CN_stem < 5
            CN_stem = 5;
        end
        CN_grain=CN_stem;
    end
    
    CNratio_leaf(i)=CN_leaf;
    CNratio_stem(i)=CN_stem;
    CNratio_grain(i)=CN_grain;
    CNratio_root(i)=CN_root;
    
    CR_ratio2(i)=CR_ratio1;
end

%%-----------------------------------------------------------------------%%
% 9. Calculate Applying CN ratio! Each component of litter has different
% C/N ratio, but they are introduced in the same soil pool!. Thus, it is
% necessary to compute the applying CN ratio! Especially when harvested!
CNratio_apply=CNratio_leaf;
CNratio_apply(HD_over_year)=(TBMLaFromLeaf(HD_over_year).*CNratio_leaf(HD_over_year)+...
    TBMLaFromStem(HD_over_year).*CNratio_stem(HD_over_year)+...
    TBMLaFromGrain(HD_over_year).*CNratio_grain(HD_over_year))...
    ./TBMLaOnlyLitterBack(HD_over_year);
NaNind=isnan(CNratio_apply);
if sum(sum(NaNind)) > 0
    CNratio_apply(NaNind)=CNratio_leaf(NaNind);
end

%%-----------------------------------------------------------------------%%
% 10. Weighted average CR_ratio and CN ratio of litterfall and root death.
% This is only used when computing kl, kh, and kd
RootDeathWhenHarvest=zeros(1,size(leaf_carbon_day_point,2));
RootDeathWhenHarvest(HD_over_year)=RootDeathWhenH;
RootdeathC_FOR_CRratio=TBMLaNoLitterBack.*CR_ratio2;
Rootdeath_For_Ave=RootdeathC_FOR_CRratio;
Rootdeath_For_Ave=Rootdeath_For_Ave+RootDeathWhenHarvest;

WeithedAveCNabove=sum(TBMLa.*CNratio_apply)./sum(TBMLa);
WeithedAveCNbelow=sum(Rootdeath_For_Ave.*CNratio_root)./sum(Rootdeath_For_Ave);
WeithedAveCNcrRatio=sum(RootdeathC_FOR_CRratio.*CR_ratio2)./sum(RootdeathC_FOR_CRratio);





















%%-----------------------------------------------------------------------%%
% Previous code






% %%%%%%%
% BD_j_day=0;
% for i=1:size(leaf_carbon_day_point,2)-1
%     % Compute Year
%     for ii=1:1:how_many_year
%         BD_detect_first_day_of_year(1)=1;
%         BD_detect_first_day_of_year(ii+1)=1+365*ii;
%     end
%     BD_detect_first_day_of_year(end)=[];
%     for ii=1:how_many_year
%         if i == BD_detect_first_day_of_year(ii)
%             BD_Check_year=BD_Check_year+1;
%             % Check How many years has been run!
%             BD_Year=start_year+BD_Check_year-1
%             BD_j_day=0;
%         end
%     end
%     % Each Year Day
%     BD_j_day=BD_j_day+1;
%     
%     % Start to compute biomass drop
%     LCi = (leaf_carbon_day_point(i))*1000; % Current step kg/m2 -> g/m2
%     %LCi = (leaf_carbon_day_point(i)+stem_carbon_day_point(i)+grain_carbon_day_point(i))*1000; % Current step kg/m2 -> g/m2
%     LCf = leaf_carbon_day_point(i+1)*1000; % Future step kg/m2 -> g/m2
%     %LCf = (leaf_carbon_day_point(i+1)+stem_carbon_day_point(i+1)+grain_carbon_day_point(i+1))*1000; % Future step kg/m2 -> g/m2
%     
%     % Compute normal turnover and phenologic changes in LAI units
%     % 60sec*60min*24hours*365days=31536000sec in a year
%     % 60sec*60min*24hours=86400sec in a day = dtime
%     if choose_crop==2
%         for j=1:size(soy_year,2)
%             if BD_Year == soy_year(j)
%                 leafspan=leafspan_Soy;
%                 %shadingF=shadingF_Soy;
%                 HavestResidueRate=HavestResidueRate_Soy;
%             end
%         end
%         for j=1:size(corn_year,2)
%             if BD_Year == corn_year(j)
%                 leafspan=leafspan_Corn;
%                 %shadingF=shadingF_Corn;
%                 HavestResidueRate=HavestResidueRate_Corn;
%             end
%         end
%     end
%     if choose_crop==3
%         if BD_Check_year == 1
%             HavestResidueRateKeep=HavestResidueRate;
%         end
%         if BD_Year == 2009
%             HavestResidueRate=HavestResidueRate2009;
%         elseif BD_Year == 2010
%             HavestResidueRate=HavestResidueRate2010;
%         else
%             HavestResidueRate=HavestResidueRateKeep;
%         end
% 
%     end
%     if choose_crop==4
%         if BD_Check_year == 1
%             HavestResidueRateKeep=HavestResidueRate;
%         end
%         if BD_Year == 2009
%             HavestResidueRate=HavestResidueRate2009;
%         else
%             HavestResidueRate=HavestResidueRateKeep;
%         end
%     end
%     
%     % 4. Abovegound biomass back to soil when harvests [%]
%     LCnt = LCi/leafspan*1;
%     
%     LCph = LCf - LCi + LCnt;
%     % Compute LAI turnover biomass
%     if LCph >= 0
%         LCTt = LCnt;
%     elseif LCph < 0
%         LCTt = (abs(LCph) + LCnt); % 
%     end
%     BMTt(i) = LCTt;
%     
%     TBMLa(i+1) = BMTt(i);
%     TBMLa(1)=BMTt(1);
% 
%     TBMLaCheck(i+1) = BMTt(i);
%     TBMLaCheck(1)=BMTt(1);
%         
%     % When harvested, put the same values as initial one
%     for j = 1:size(HD_over_year,2)
%         if i  == HD_over_year(j)
%             TBMLa(i) = TBMLa(1);
%             
%             % For checking
%             TBMLaCheck(i) = TBMLa(1);
%             TBMLaCheck(i) = Tot_AboveC(i).*HavestResidueRate;
%             
%             % TBMLa only Litter back when harvest
%             TBMLaOnlyLitterBack(i)=Tot_AboveC(i).*HavestResidueRate;
%             TBMLaFromLeaf(i) = leaf_carbon_day_point(i).*1000.*HavestResidueRate;
%             TBMLaFromStem(i) = stem_carbon_day_point(i).*1000.*HavestResidueRate;
%             TBMLaFromGrain(i) = 0;
%             if choose_crop == 3
%                 TBMLaFromGrain(i) = grain_carbon_day_point(i).*1000.*HavestResidueRate;
%             end
%             
%         end
%     end
% end
% TBMLaNoLitterBack=TBMLa;
% 
% 
% % For Plotting cumulative TBMLa
% for i=1:size(BD_detect_first_day_of_year,2)
%     CumulTBMLaCheck=cumsum(TBMLaCheck(365*(i-1)+1:365*i));
%     CumulTBMLaEachYearCheck(365*(i-1)+1:365*i)=CumulTBMLaCheck;
% end
% 
% % This is the Juan's Eqn
% TBMLa = TBMLaCheck;