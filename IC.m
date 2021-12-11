% Set initial conditions.

function [u1,x] = IC(i1,i2,N,dz)
 pi = 2*acos(0);
 lambda = 1;
%     for i = i1:i2
%         x(i) = dx*(i-i1);
%         u1(i) = sin(2.*pi/N*x(i));
%     end
% x = zeros(length(dz),1);
    for i = i1+1:i2
%         x(i) = x(i-1) + dz(i);
%         if (x(i)>=20) && (x(i)<=70)
%             u1(i)   = exp(-0.01*(x(i)-45)^2);
%         end
%         u1(i)   = exp(-Kx*pi^2*t)*sin(pi*x(i));  

%         if (x(i)>=10) && (x(i)<=30)
%             u1(i)=1.0*exp(-0.05*(x(i)-20)^2); %gaussian
%         else u1(i) = 0;
%         end
%         u1(i) = sin(2.*pi/N*x(i));
%         u1(i)=lambda*exp(-lambda/20*dz(i)^2);
        u1(i)=lambda*exp(-lambda*dz(i)^2);
    end
    u1(i2+1) = u1(i2);
%     figure()
%         plot(u1,dz)
%         xlabel('Concentration');
%         ylabel('depth');
%         title('Initial Condition');
%         set(gca,'YDir','reverse');
end