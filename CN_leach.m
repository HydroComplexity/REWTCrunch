% modified his code
%function [LCH_CC, TLCH_CC] = CN_leach(PARAMS, CC, aa, qq, sm, dz_mm)
function [LCH_CC, TLCH_CC, flow_CC_m2, LCH_CC_Top_Up, LCH_CC_Bottom_Down, LCH_CC_Top_Down, LCH_CC_Bottom_Up]...
    = CN_leach(PARAMS, CC, aa, qq, sm, dz_mm)

    % CC Concentration variable [gr/m3]
    % aa solubility of element in study
    nl_soil=PARAMS.nl_soil;%       nl_soil = # soil layers
    leach_on = 1;
%     leach_on = PARAMS.CN.leach_on;%       leach_on = Leaching ON or OF [1 0]
    CC(nl_soil+1) = 0; %for upward flux calc. SUSANA EDIT
    dz_mm(nl_soil+1) = dz_mm(nl_soil); %for upward flux calc. SUSANA EDIT
    sm(nl_soil+1) = sm(nl_soil);
    % Declare matrix
    LCH_CC = zeros (nl_soil,1);
    
    LCH_CC_Top_Up      = zeros (nl_soil,1);
    LCH_CC_Top_Down    = zeros (nl_soil,1);
    LCH_CC_Bottom_Up   = zeros (nl_soil,1);
    LCH_CC_Bottom_Down = zeros (nl_soil,1);
    
    AA = zeros (nl_soil,1);
    BB = zeros (nl_soil,1);
%     LCH_CC_Age_Top_Up=zeros (size(CC_Age));
%     LCH_CC_Age_Top_Down=zeros (size(CC_Age));
%     LCH_CC_Age_Bottom_Up=zeros (size(CC_Age));
%     LCH_CC_Age_Bottom_Down=zeros (size(CC_Age));
%     CC_Age_m2(:,:,1)=CC_Age(:,:,1).*repmat(dz_mm,1,size(CC_Age,2))./1000;
%     CC_Age_m2(:,:,2)=CC_Age(:,:,2).*repmat(dz_mm,1,size(CC_Age,2))./1000;
%     CC_Age_m2(:,:,3)=CC_Age(:,:,3).*repmat(dz_mm,1,size(CC_Age,2))./1000;
%     Tolerance_Age=10^-10;
    % dk
    %for ii = 1:nl_soil  
    %   if sm(ii) == 0
    %      sm(ii)=9.9999*10^(-7);
    %   end
    %end
    
    
% below is Juan' code    
    for ii = 1:nl_soil    
        %   positive value = loss (subtracted from budget)
        %   NOTE: first H2O flux value is from infiltration into soil column
            LCH_CC(ii) = 0;
            if (leach_on)   
                if (ii>1)   % transport through top layer interface
                    if (qq(ii)>0)  % downward transport into cell -
                        flux=(aa/sm(ii-1)*CC(ii-1)*qq(ii));
                        if flux > aa*CC(ii-1)*(dz_mm(ii-1)/1000)
                            flux = aa*CC(ii-1)*(dz_mm(ii-1)/1000);
                        end
                        
                        %LCH_CC(ii) = LCH_CC(ii) - (aa/sm(ii-1)*CC(ii-1)*qq(ii));% [gr/m^2/d] 
                        %LCH_CC_Top_Down(ii) = - (aa/sm(ii-1)*CC(ii-1)*qq(ii));
%                         LCH_CC(ii) = LCH_CC(ii) - (aa/sm(ii-1)*CC(ii-1)*qq(ii));% [gr/m^2/d] %susana
%                     commented
                        LCH_CC_Top_Down(ii) = - flux; %- (aa/sm(ii-1)*CC(ii-1)*qq(ii));    %susana edit: replaced "flux"                    
                        
                    else           % upward transport out of cell -
                        flux=(aa/sm(ii)*CC(ii)*qq(ii));
                        if flux < -aa*CC(ii)*(dz_mm(ii)/1000)
                            flux = -aa*CC(ii)*(dz_mm(ii)/1000); %susana edit: (-) sign added
                        end
                        
                        
                        %LCH_CC(ii) = LCH_CC(ii) - (aa/sm(ii)*CC(ii)*qq(ii));% [gr/m^2/d]  
                        %LCH_CC_Top_Up(ii)=- (aa/sm(ii)*CC(ii)*qq(ii));
%                         LCH_CC(ii) = LCH_CC(ii) - flux;% [gr/m^2/d]  %susana
%                     commented
                        LCH_CC_Top_Up(ii)=- flux;                        
                        
                    end
                end
                
            %    % transport through bottom layer interface
                if (qq(ii+1)>0)  % downward transport out of cell +
                    flux=(aa/sm(ii)*CC(ii)*qq(ii+1));
                    if flux > aa*CC(ii)*(dz_mm(ii)/1000)
                       flux = aa*CC(ii)*(dz_mm(ii)/1000);
                    end
                    %LCH_CC(ii) = LCH_CC(ii) + (aa/sm(ii)*CC(ii)*qq(ii+1));% [gr/m^2/d]
                    %LCH_CC_Bottom_Down(ii) = (aa/sm(ii)*CC(ii)*qq(ii+1));% [gr/m^2/d]
%                     LCH_CC(ii) = LCH_CC(ii) + flux;% [gr/m^2/d] %susana
%                     commented
                    LCH_CC_Bottom_Down(ii) = flux;% [gr/m^2/d]                    

                else           % upward transport into cell +
                    flux=(aa/sm(ii+1)*CC(ii+1)*qq(ii+1));
                    if flux < -aa*CC(ii+1)*(dz_mm(ii+1)/1000)
                       flux = -aa*CC(ii+1)*(dz_mm(ii+1)/1000); %susana edit: added (-)
                    end
                    
                    %LCH_CC(ii) = LCH_CC(ii) + (aa/sm(ii+1)*CC(ii+1)*qq(ii+1));% [gr/m^2/d]
                    %LCH_CC_Bottom_Up(ii) = (aa/sm(ii+1)*CC(ii+1)*qq(ii+1));% [gr/m^2/d]
%                     LCH_CC(ii) = LCH_CC(ii) + flux;% [gr/m^2/d] %susana
%                     commented
                    LCH_CC_Bottom_Up(ii) = flux;% [gr/m^2/d]                    

                end
                if (ii==nl_soil)
                    flux=(aa/sm(ii)*CC(ii)*qq(ii+1));
                    if flux > aa*CC(ii)*(dz_mm(ii)/1000)
                       flux = aa*CC(ii)*(dz_mm(ii)/1000);
                    end
                    %TLCH_CC = (aa/sm(ii)*CC(ii)*qq(ii+1));% [gr/m^2/d]
                    TLCH_CC = flux;% [gr/m^2/d]
                end   
            end                        
    end
    
    %LCH_CC_Bottom_Up+LCH_CC_Bottom_Down+LCH_CC_Top_Up+LCH_CC_Top_Down
    % Check the error
%     mberrorLCH=sum(sum((sum(sum(LCH_CC_Age_Top_Up(:,:,:),2),3)+sum(sum(LCH_CC_Age_Top_Down(:,:,:),2),3)...
%         +sum(sum(LCH_CC_Age_Bottom_Up(:,:,:),2),3)+sum(sum(LCH_CC_Age_Bottom_Down(:,:,:),2),3))...
%         -(LCH_CC)));
    
    % Dongkook Error Correct!
    TLCH_CC = TLCH_CC(:);
    %LCH_CC = LCH_CC(:);   
    LCH_CC=LCH_CC_Bottom_Up+LCH_CC_Bottom_Down+LCH_CC_Top_Up+LCH_CC_Top_Down;
    
    % When qq i) top layer goes up ii) bottom layer goes down. Negative
    % occurs. Fix
    %CC_track_m2=CC_track.*(dz_mm./1000);
    %indNeg=(CC_track_m2-LCH_CC)<0;
    %AA=LCH_CC_Bottom_Down;
    %BB=LCH_CC_Top_Up;
    %LCH_CC_Bottom_Down(indNeg)=(AA(indNeg)/(AA(indNeg)+BB(indNeg)))*CC_track_m2(indNeg);
    %LCH_CC_Top_Up(indNeg)=(BB(indNeg)/(AA(indNeg)+BB(indNeg)))*CC_track_m2(indNeg);
      
    % susana edit begin
%     LCH_CC_Top_Down(2:end)=-LCH_CC_Bottom_Down(1:end-1);
%     LCH_CC_Bottom_Up(1:end-1)=-LCH_CC_Top_Up(2:end);
    
%     LCH_CC=LCH_CC_Bottom_Up+LCH_CC_Bottom_Down+LCH_CC_Top_Up+LCH_CC_Top_Down;
    % susana edit end
    
    % Check Point
%     CCnew = CC.*(dz_mm./1000) - (LCH_CC); % g/m2
%     ind = CCnew < 0;
    %CCnew(ind) = 0;

    
    
    
    %LCH_CC_Bottom_Up(1:end-1,:)=-(LCH_CC_Top_Up(2:end,:));
    %LCH_CC_Top_Down(2:end,:)=-(LCH_CC_Bottom_Down(1:end-1,:));
    
%     mberrorLCH=sum((LCH_CC_Bottom_Up+LCH_CC_Bottom_Down+LCH_CC_Top_Up+LCH_CC_Top_Down)-LCH_CC);
    
    
    
    
    
    
    
    
    
    
    
    % dk to calclate only flow from the bottom of each layer.
    for ii = 2:nl_soil
        if (qq(ii)>0) % g/m2
            flow_CC_m2(ii-1)=(aa/sm(ii-1)*CC(ii-1)*qq(ii));
        else
            flow_CC_m2(ii-1)=(aa/sm(ii)*CC(ii)*qq(ii));
        end
    end
    flow_CC_m2=flow_CC_m2(:);
            
%            if (leach_on)
%                if (ii>1)   % transport through top layer interface
%                    if (qq(ii)>0)  % downward transport into cell
%                        LCH_CC(ii) = LCH_CC(ii) - (aa/sm(ii-1)*CC(ii-1)*qq(ii));% [gr/m^2/d]
%                    else           % upward transport out of cell
%                        LCH_CC(ii) = LCH_CC(ii) - (aa/sm(ii)*CC(ii)*qq(ii));% [gr/m^2/d]
%                    end
%                end
%                % transport through bottom layer interface
%
%                if (ii==nl_soil && qq(ii+1)<0)
%                        LCH_CC(ii) = 0;% [gr/m^2/d]
%                        TLCH_CC = 0;% [gr/m^2/d]                
%                else
%                    if (qq(ii+1)>0)  % downward transport out of cell
%                        LCH_CC(ii) = LCH_CC(ii) + (aa/sm(ii)*CC(ii)*qq(ii+1));% [gr/m^2/d]
%                    else           % upward transport into cell
%                        LCH_CC(ii) = LCH_CC(ii) + (aa/sm(ii+1)*CC(ii+1)*qq(ii+1));% [gr/m^2/d]
%                    end
%                    TLCH_CC = (aa/sm(ii)*CC(ii)*qq(ii+1));% [gr/m^2/d]
%                end   
%            end                        
%    end

% until here

% dk's problem in mberrorN
%    for ii = 1:nl_soil    
%        %   positive value = loss (subtracted from budget)
%        %   NOTE: first H2O flux value is from infiltration into soil column
%           LCH_CC(ii) = 0;
%           if (leach_on)
%                if (ii>1)   % transport through top layer interface
%                    if (qq(ii)>0)  % downward transport into cell
%                        LCH_CC(ii) = LCH_CC(ii) - (CC(ii-1)*(qq(ii)/(qq(ii)+(sfield(ii-1))*dz_mm(ii-1)/1000)));% [gr/m^2/d]
%                    else           % upward transport out of cell
%                        
%                        LCH_CC(ii) = LCH_CC(ii) - (CC(ii)*(qq(ii)/(qq(ii)+(sfield(ii))*dz_mm(ii)/1000)));% [gr/m^2/d]
%                    end
%                end
%                % transport through bottom layer interface
%                if (qq(ii+1)>0)  % downward transport out of cell
%                    LCH_CC(ii) = LCH_CC(ii) + (CC(ii)*(qq(ii+1)/(qq(ii+1)+(sfield(ii))*dz_mm(ii)/1000)));% [gr/m^2/d]
%                else           % upward transport into cell
%                    LCH_CC(ii) = LCH_CC(ii) + (CC(ii+1)*(qq(ii+1)/(qq(ii+1)+(sfield(ii+1))*dz_mm(ii+1)/1000)));% [gr/m^2/d]
%                end
%                if (ii==nl_soil)
%                    TLCH_CC = (CC(ii)*(qq(ii+1)/(qq(ii+1)+(sfield(ii))*dz_mm(ii)/1000)));% [gr/m^2/d]
%                end   
%            end                        
%    end



    