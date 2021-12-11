function [ puptot ] = CN_reuptake( cexremain, theta, wuptake_all_store, dz,i1,i2,n,extype)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


rsd = 0.1; %rescaled diffusion coefficient;
nu = 0.6; %typically between 0.5 and 0.8

if extype == 1
    % aa = solubility of glucose;
%     aa = 909; %mg/mL at 25C (online, incl Wikipedia:( )
%     aa = aa*10^(-3); %g/mL
    aa = 1; % removing this 7/11/19 because what I needed was fraction dissolved (lambda), not solubility.
    lambda = 0.7; %NEED TO FIND A REFERENCE FOR THIS
    for i = i1:i2
        pup_w(i)    = aa*wuptake_all_store(i)./dz(i)*1000*lambda*cexremain(i)/theta(i);
        ku(i)       = 1./(theta(i).*n(i).*dz(i))*rsd.*theta(i)^nu;
        pup_d(i)    = ku(i)*lambda*cexremain(i);
        puptot(i)   = pup_w(i) + pup_d(i);
    end
elseif extype == 2
    for i = i1:i2        
        puptot(i)   = 0;
    end
end




end

