function [ u2amv, u2advup, u2advdown ] = CN_exadvection( u1, i1, i2, v, dz, theta, extype )

u2advup = zeros(i2,1);
u2advdown = zeros(i2,1);

if extype == 1
    % aa = solubility of glucose;
    aa = 909; %mg/mL at 25C (online, incl Wikipedia:( )
    aa = aa*10^(-3); %g/mL
elseif extype == 2
    % aa = solubility of flavonoids;
    aa = 0.00215; %g/L at 25C (sciencedirect, reported for querticin)
    aa = aa*10^(-3); %g/mL
end

for ii = i1:i2    
        %   positive value = loss (subtracted from budget)
        %   NOTE: first H2O flux value is from infiltration into soil column
            u2amv(ii) = 0;
                if (ii>1)   % transport through top layer interface
                    if (v(ii)>0)  % downward transport into cell
                        u2amv(ii) = u2amv(ii) - (aa/theta(ii-1)*u1(ii-1)*v(ii));% [gr/m^2/d]
                        if isinf(u2amv(ii))
                            u2amv(ii)=0;
                        end
                    else           % upward transport out of cell
                        u2amv(ii) = u2amv(ii) - (aa/theta(ii)*u1(ii)*v(ii));% [gr/m^2/d]
                        u2advup(ii) = - (aa/theta(ii)*u1(ii)*v(ii));
                        if isinf(u2amv(ii))
                            u2amv(ii)=0;
                        end
                    end
                end
                % transport through bottom layer interface
                if (v(ii+1)>0)  % downward transport out of cell
                    u2amv(ii) = u2amv(ii) + (aa/theta(ii)*u1(ii)*v(ii+1));% [gr/m^2/d]
                    u2advdown(ii) =  (aa/theta(ii)*u1(ii)*v(ii+1));
                    if isinf(u2amv(ii))
                        u2amv(ii)=0;
                    end
                else           % upward transport into cell
                    u2amv(ii) = u2amv(ii) + (aa/theta(ii+1)*u1(ii+1)*v(ii+1));% [gr/m^2/d]
                    if isinf(u2amv(ii))
                        u2amv(ii)=0;
                    end
                end
                if (ii==i2)
                    u2anomv = (aa/theta(ii)*u1(ii)*v(ii+1));% [gr/m^2/d]
                end   

%     u2a(ii+1) = - v(ii)./dz(ii)*(u1(ii+1)-u1(ii));
end
    u2anomv = u2anomv(:);
    u2amv = u2amv(:);
    
%     CCnew = CC - LCH_CC;
%     ind = CCnew < 0;
%     CCnew(ind) = 0;

end

