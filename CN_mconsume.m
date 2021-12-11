function [ mup,Cbglu,cexremain ] = CN_mconsume( i1,i2,Cb,cex,extype)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

aflav = 0.85; % fraction of flavonoid concentration affecting microbial activity.
% starve = zeros(size(cex));
    if extype == 1
        rg = 0.17*10^(-6); %g glucose/g microbes/hour 
        rg = rg*24; %g glucose/g microbes/day
        for i = i1:i2
            gmuppot(i) = rg*Cb(i); % g glucose consumed/day
            if cex(i) > 0
                if cex(i) <= gmuppot(i) %if the cex is less than what all of the microbes need to consume, some will die
                    mup(i) = cex(i);
                    starve(i) = gmuppot(i) - mup(i);
                    cexremain(i) = 0;  
                else 
                    mup(i) = gmuppot(i);       
                    cexremain(i) = gmuppot(i)-cex(i);
                    starve(i) = 0;
                end
           
                Cbglu = 0.39.*(mup-starve).*0.4; %gmup*0.4=C consumed from glucose; 0.39xthat is additional C biomass added from consumption (Behera and Wagner 1974)
                Cbglu = Cbglu';
            %     Cbglu = [Cbglu; 0]; %adding for last flux outside; other variables have 13 layers because they're fluxes
            else 
                Cbglu = zeros(size(cex))';
                mup = zeros(size(cex));
                cexremain = zeros(size(cex));
            end
        end
    elseif extype == 2
        for i = i1:i2
            mup(i) = aflav*cex(i);
            cexremain(i) = (1-aflav)*cex(i);
            Cbglu = 0;
        end
    end

    
end

