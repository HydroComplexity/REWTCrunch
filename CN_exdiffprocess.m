function [ u2d, uup_sub, udown_sub ] = CN_exdiffprocess( u1,i1,i2,dz,extype,cex)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
udown = zeros(1,i2);
uup = udown;
u2d = udown;
udown_sub = udown;
uup_sub = udown;

if extype == 1
    D = 0.665*10^-5; % Diffusion coeff. for glucose [cm2/s]
    D = D/(100^2)*24*60*60; %changes from cm2/s to m2/day
elseif extype == 2
    D = 1.89*10^(-7); % Diffusion coeff. for mixture of exudates (assuming includes flav) [cm2/s]
    D = D/(100^2)*24*60*60; %changes from cm2/s to m2/day;
end

for ii=i1:i2-1 
%  		u2d(ii+1) = D./dz(ii)^2*(u1(ii+2) - 2.*u1(ii+1) + u1(ii));	
    if u1(ii) > u1(ii+1)
        %udown(ii+1) = -D*(u1(ii)-u1(ii+1))/((dz(ii+1)+dz(ii))/2); %[g/m2/d]
        udown_sub(ii) = cex(ii)*D*(u1(ii)-u1(ii+1))/((dz(ii+1)+dz(ii))/2); %[g/m2/d]            
    else
        %uup(ii) = -D*(u1(ii+1)-u1(ii))/((dz(ii+1)+dz(ii))/2); %[g/m2/d]
        uup_sub(ii+1) = cex(ii+1)*D*(u1(ii+1)-u1(ii))/((dz(ii+1)+dz(ii))/2); %[g/m2/d]
    end
end
    
%     for i = i1:i2-1
%         if i == i1
%             u2d(i) = -udown(i) + uup(i+1);
%         else
%             u2d(i) = udown(i-1) - udown(i) - uup(i) + uup(i+1);
%         end       
%     end
    
   
    
    
    
    
%     u2d(i2) = u2d(i2-1); %for last layer

%     u2d = u2d';

end

