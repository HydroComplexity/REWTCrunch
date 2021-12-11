function [pup_new,mup_new,u2advup_new,u2advdown_new,u2advup_add_new,u2advdown_add_new,...
    uup_new,udown_new,uup_add_new,udown_add_new] = FixNegative(u1,i1,i2,pup,mup,u2advup_m2,u2advdown_m2,uup_m2,udown_m2,dz)
% Total loss from layer
%gpup and gmup are positive
%adv positive outward
pup_m2 = pup'.*dz;
mup_m2 = mup'.*dz;
u1_m2 = u1(i1:i2)'.*dz;

u2advup_m2_add_new = zeros(i2,1);
u2advdown_m2_add_new = zeros(i2,1) ;
uup_m2_add_new =  zeros(i2,1);
udown_m2_add_new =  zeros(i2,1);
exleach_add_new =  zeros(i2,1);

       
for i = i1:i2
%     Totalloss(i) = pup_m2(i) + mup_m2(i) + u2advup_m2(i) + u2advdown_m2(i) + uup_m2(i) + udown_m2(i) + exleach(i);
   Totalloss(i) = pup_m2(i) + mup_m2(i) + u2advup_m2(i) + u2advdown_m2(i) + uup_m2(i) + udown_m2(i);
   if Totalloss(i) > u1_m2(i)
       gpup_m2_new(i) = u1_m2(i).*(pup_m2(i)./Totalloss(i));
       gmup_m2_new(i) = u1_m2(i).*(mup_m2(i)./Totalloss(i));
       u2advup_m2_new(i) = u1_m2(i).*(u2advup_m2(i)./Totalloss(i));
       u2advdown_m2_new(i) = u1_m2(i).*(u2advdown_m2(i)./Totalloss(i));
       uup_m2_new(i) = u1_m2(i).*(uup_m2(i)./Totalloss(i));
       udown_m2_new(i) = u1_m2(i).*(udown_m2(i)./Totalloss(i));
%        exleach_new(i) = u1_m2(i).*(exleach(i)./Totalloss(i));
    else
       gpup_m2_new(i) = pup_m2(i);
       gmup_m2_new(i) = mup_m2(i);
       u2advup_m2_new(i) = u2advup_m2(i);
       u2advdown_m2_new(i) = u2advdown_m2(i);
       uup_m2_new(i) = uup_m2(i);
       udown_m2_new(i) = udown_m2(i);
%        exleach_new(i) = exleach(i);
        
    end
end

% Adv and disp adding to the layers
u2advup_m2_add_new(1:end-1) = -u2advup_m2_new(2:end);
u2advdown_m2_add_new(2:end) = -u2advdown_m2_new(1:end-1);
uup_m2_add_new(1:end-1) =  -uup_m2_new(2:end);
udown_m2_add_new(2:end) =  -udown_m2_new(1:end-1);
% exleach_add_new(2:end) =  -exleach_new(1:end-1);

% Conversion
pup_new = pup_m2./dz;
mup_new = mup_m2./dz;
u2advup_new = u2advup_m2_new'./dz;
u2advdown_new = u2advdown_m2_new'./dz;
uup_new = uup_m2_new'./dz;
udown_new = udown_m2_new'./dz;
u2advup_add_new = u2advup_m2_add_new./dz;
u2advdown_add_new = u2advdown_m2_add_new./dz;
uup_add_new =  uup_m2_add_new./dz;
udown_add_new = udown_m2_add_new./dz;
% exleach_new = exleach_new'./dz;
end
