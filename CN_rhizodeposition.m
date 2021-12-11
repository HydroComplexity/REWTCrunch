function [ cex, Bp, Br, tn ] = CN_rhizodeposition(doy,day,iscorn,rootfr,i1,i2,extype,LAIday,dz,BulkDensity,Amm,Nit,te,th,VARIABLES,Br_print)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%      dz = dz./1000; %back to mm
  Brprev = VARIABLES.SOIL.Br;
% Brprev = Br_print(:,day-1);
  
   if iscorn == 1 
        % Get total root biomass for plant
        if doy(day) < te-1
            Bp = 0;
            %rp = 0.8; % R:S (Bonifas 2005)
            rp = -0.0079*(doy(day)-te) + 0.8266;%(Bonifas 2005, Fig1)
            Br = rp*Bp; % root biomass g/m2          

        elseif doy(day) > te-1 &&  doy(day) < th+1
            tn = doy(day) - te;
            Bp = 426.21*LAIday - 23.245; %Gilabert, et al. 1996
            if Bp < 0
                Bp = 0;
            end
            %Bp = LAIday*13.48*(6*10^(-7)*tn^4 - 0.0001*tn^3 + 0.0095*tn^2 - 0.2567*tn + 2.8); % g/m2
%             Bp = Bp*0.4127; % converts g/m2 to gC/m2 (Latshaw 1924)
            rp = -0.0079*(doy(day)-te) + 0.8266;%(Bonifas 2005, Fig1)
            if rp < 0.28
                rp = 0.27; %smallest value at maturity (Bonifas)
            end
            %rp = 0.8; % R:S (Bonifas 2005)
            Br = rp*Bp; % root biomass g/m2                
        elseif doy(day) == th+1 %harvest
            Bp = 345/0.4127; %avg g/m2 (Anderson-Teixeira 2013)/(Latshaw 1924)
            rp = -0.0079*(doy(day)-te) + 0.8266;%(Bonifas 2005, Fig1)
%             rp = 0.27; % R:S (Bonifas 2005)
            if rp < 0.28
                rp = 0.27; %smallest value at maturity (Bonifas)
            end
            Br = rp*Bp; % root biomass g/m2   
        else
            Bp = 0; %post harvest
            Br = Brprev - 0.21*Brprev; % previous root biomass - death rate    
            if Br < 0
                Br = 0;
            end
        end 

        if extype == 1 % Glucose  
%             rex = rex(1);
            rex = 4.82*10^-3; %g glucose/ g dry root /day 
%            rex = 0;
            for i = i1:i2        
                rho_r(i) = rootfr(i).*Br; %     rootfr is in percent
                cex(i) = rex.*rho_r(i);
            end
        elseif extype == 2    % Flavonoid    
            rex = 50*10^-9; % mol/h/tip ; can be 50-225
            rex = rex*24; % mol flavonoid/day/tips these values are too
%             high for the BNI conversion...units must be off
%               rex = rex(2);
%               rex = 2.5*10^-9; %mol/h/tip (Kidd 2001)
%               rex = rex*24;
%             rex = 0;
            ammtonit = 0.06; % minimum ratio of ammonium to nitrogen before flav kicks in (based on Subbarao 2013)
            for i = i1:i2
                if Br > 0
                    rtips(i) = Br*rootfr(i); % root tips per layer  (est. 1000 root tips per plant; Gerson Meneghetti Sarzi Sartori 2016)
                else
                    rtips(i) = 0;
                end
                flavthresh(i) = (Amm(i)-Nit(i))/Nit(i);
                if flavthresh(i) <= ammtonit
                    cex(i) = rex*rtips(i); %mol flav/day
                else
                    cex(i) = 0;
                end
            end
        end
   else %SOYBEAN
        % Get total root biomass for plant
        % max Br for soybean is about 250g;
       
        if doy(day) < te-1
            Bp = 0;
            rp = 0.23; % R:S (Rogers 1992) under ambient CO2            
            Br = rp*Bp; % root biomass g/m2    
        elseif doy(day) > te-1 &&  doy(day) < th+1
            tn = doy(day) - te;
            Bp = 31.11*LAIday^2 - 32.95*LAIday; % g/m2%     
            if Bp < 0
                Bp = 0;
            end
%             Bp = Bp*0.3520; %conversion to gC/m2 (Srivastava 2006)
            rp = 0.23; % R:S (Rogers 1992) under ambient CO2 
            Br = rp*Bp; % root biomass g/m2    
        else
            Bp = 0;
            Br = Brprev - 0.35*Brprev; % previous root biomass - death rate    
            if Br < 0
                Br = 0;
            end
        end 

        if extype == 1 % Glucose 
            if Br > 0
%                 rex = rex(1);
                rex = 42*10^(-6); %g glucose/ root network /15 min (Timotiwu and Sakurai 2002) NEED TO FIX THIS
                rex = rex/15*60*24; %g glucose/root network/ day
                rex = Br/250*rex;
%                 rex = 0;
            else
                rex = 0;
            end
            for i = i1:i2        
                cex(i) = rex.*rootfr(i);
            end
        elseif extype == 2    % Flavonoid    
%             rex = 1000*10^(-9); % mol/g root (Graham 1991)
%             rex = rex(2);
            rex = 1200*10^(-12); % mol/root (Peuppke 1998)
%             rex = rex*1000; % mmol/root
%             rex = 0.2*10^(-6); % g/ g dry soil / 7min (Guo, et al 2011)
%             rex = rex/7*60*24; % g/ g dry soil/ day
%             rex = rex/270.24/1000; % mmol/ g dry soil/day (MW_genistein = 270.24 g/mol (Pubchem))
%             rex = 0;
            
            ammtonit = 0.06; % minimum ratio of ammonium to nitrogen before flav kicks in
            
            for i = i1:i2
                flavthresh(i) = (Amm(i)-Nit(i))/Nit(i);
                if flavthresh(i) <= ammtonit && Br > 0
                    if Br > 0
                        rtips(i) = floor(Br*10*rootfr(i)); % root tips per layer  (est. 3700 root tips per plant; Wu and Guo 2014)
                    else
                        rtips(i) = 0;
                    end
%                     cex(i) = rex.*dz(i).*BulkDensity(i)*rootfr(i); % mmol/day
%                     cex(i) = rex.*Br*rootfr(i);
                    cex(i) = rex*rtips(i); %mol flav/day
                else
                    cex(i) = 0;
                end
            end
        end
   end
   


end

