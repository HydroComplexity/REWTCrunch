function [gmup,Cbglu,Gtot,Ftot,fflav,cexf,cexg, Bp, Br,VARIABLES] = ...
    CN_exudate(doy,day,iscorn,VARIABLES,VERTSTRUC,qq,layeruptake_all,PARAMS,smfull,rootfr,BulkDensity,te,th)

% extype          = [1,2];
% extypename      = {'Glucose','Flavonoids'};
dz              = VERTSTRUC.dzsmm./1000; % m
i1              = 1;
i2              = PARAMS.nl_soil;
n               = VARIABLES.timestep;
qlayer          = qq;
v               = qlayer; % [m/d]
v(i1)           = 0; % first layer ET
%v(i2+1)         = v(i2); % last layer flux boundary
Cb              = VARIABLES.Cb;
theta           = smfull;
n               = PARAMS.porosity;
ndz             = length(dz); % number of soil layers
N               = ndz;
layeruptake     = layeruptake_all;
hk(i1:i2)       = VARIABLES.SOIL.hk;
C_amm           = VARIABLES.Amm_mo;
fflav           = 0;
lambda          = 0.8; 
Amm             = VARIABLES.Amm;
Nit             = VARIABLES.Nit;

if iscorn == 1
%     yrLAI           = load('lai1yrcorn.mat'); % from MODIS filling
%     LAIday          = yrLAI.lai1yrcorn.*1.5; %temporary mult before using actual numbers from Landsat
    LAIday = cell2mat(struct2cell(load('smoothlaicorn.mat'))); %smoothed modis
    LAIday(1:te-1)=0;
    LAIday(th+1:end)=0;
else
%     yrLAI           = load('lai1yrsoy.mat'); % from MODIS filling; temp til Landsat
%     LAIday          = yrLAI.lai1yr;
    LAIday = cell2mat(struct2cell(load('smoothlaisoy.mat'))); %smoothed modis
    LAIday(1:te-1)=0;
    LAIday(th+1:end)=0;
end
    

if doy(day) == 1
    cexgprev = 0;
    VARIABLES.SOIL.cexgprev = cexgprev;
    cexfprev = 0;
    VARIABLES.SOIL.cexfprev = cexfprev;
    Brprev = 0;
    VARIABLES.SOIL.Br = Brprev;
else
    cexgprev = VARIABLES.SOIL.cexgprev;
    cexfprev = VARIABLES.SOIL.cexfprev;
end

% % Initial Conditions
    [u1] = IC(i1,i2,N,dz);


% Exudation
LAIday = LAIday(doy(day));
% LAIday = 0;
[cexg, Bp, Br] = CN_rhizodeposition(doy,day,iscorn,rootfr,i1,i2,1,LAIday,dz,BulkDensity,Amm,Nit,te,th,VARIABLES); %g glucose/day
% cexg = cexg+cexgprev;
[cexf, Bp, Br] = CN_rhizodeposition(doy,day,iscorn,rootfr,i1,i2,2,LAIday,dz,BulkDensity,Amm,Nit,te,th,VARIABLES); % for flavonoids
% cexf = cexf+cexfprev;
% Bp = zeros(size(Bp));
% Br = zeros(size(Br));
VARIABLES.SOIL.Br = Br;
% cexg = zeros(size(cexg));
% cexf = zeros(size(cexf));

% Advection
    
%     u1 = BC(i1,i2,u1);
% The unit of adv is [g/m2/timestep]
% Positive = outward
[u2ag, u2advup_m2g, u2advdown_m2g] = CN_exadvection(u1, i1, i2, v, dz, theta,1); % for glucose
[u2af, u2advup_m2f, u2advdown_m2f] = CN_exadvection(u1, i1, i2, v, dz, theta,2); % for flavonoids
% u2ag = zeros(size(u2ag)); 
% u2advup_m2g = zeros(size(u2advup_m2g));
% u2advdown_m2g = zeros(size(u2advdown_m2g));
% u2af = zeros(size(u2af)); 
% u2advup_m2f = zeros(size(u2advup_m2f));
% u2advdown_m2f = zeros(size(u2advdown_m2f));

% Diffusion  
% Diff unit: [g/m2/timestep]
% All negative number adding to the layer
% There are 2 variables: Up and Down
[u2dg, uup_m2g, udown_m2g] = CN_exdiffprocess(u1,i1,i2,dz,1,cexg); % for glucose
[u2df, uup_m2f, udown_m2f] = CN_exdiffprocess(u1,i1,i2,dz,2,cexf); % for flavonoids
% u2dg = zeros(size(u2dg)); 
% uup_m2g = zeros(size(uup_m2g));
% udown_m2g = zeros(size(udown_m2g)); 
% u2df = zeros(size(u2df)); 
% uup_m2f = zeros(size(uup_m2f));
% udown_m2f = zeros(size(udown_m2f));

% Microbial activity
[gmup, Cbglu,cexgremain] = CN_mconsume(i1,i2,Cb,cexg,1);
[fmup, Cbf,cexfremain] = CN_mconsume(i1,i2,Cb,cexf,2);
% Cbglu = zeros(size(Cbglu));
% Cbf = zeros(size(Cbf));

% BNI factor for flavonoids
[fflav] = CN_BNIfactor(i1,i2,cexf);

% Plant Uptake - assumed to be zero for flavonoids
gpup = CN_reuptake(cexgremain, theta, layeruptake, dz,i1,i2,n,1); %g/day
fpup = CN_reuptake(cexfremain, theta, layeruptake, dz,i1,i2,n,2); %g/day

% Leaching  
% [gexleach] = CN_exleach(i1,i2,theta,n,cexgremain,1,v,dz);
% [fexleach] = CN_exleach(i1,i2,theta,n,cexfremain,2,v,dz);
% 


% Fix negative solution 
[gpup_new,gmup_new,u2advupg_new,u2advdowng_new,u2advupg_add_new,u2advdowng_add_new,...
    uupg_new,udowng_new,uupg_add_new,udowng_add_new] = ...
    FixNegative(u1,i1,i2,gpup,gmup,u2advup_m2g,u2advdown_m2g,uup_m2g,udown_m2g,dz); 

[fpup_new,fmup_new,u2advupf_new,u2advdownf_new,u2advupf_add_new,u2advdownf_add_new,...
    uupf_new,udownf_new,uupf_add_new,udownf_add_new] = ...
    FixNegative(u1,i1,i2,fpup,fmup,u2advup_m2f,u2advdown_m2f,uup_m2f,udown_m2f,dz); 

% Total
Gtot = u1(i1:i2)' + cexg' - gpup_new - gmup_new ...
    - (u2advupg_new + u2advdowng_new + u2advupg_add_new + u2advdowng_add_new) ...
    - (uupg_new + udowng_new + uupg_add_new + udowng_add_new); %change in at every time step

Ftot = u1(i1:i2)' + cexf' - fpup_new - fmup_new ...
    - (u2advupf_new + u2advdownf_new + u2advupf_add_new + u2advdownf_add_new) ...
    - (uupf_new + udownf_new + uupf_add_new + udownf_add_new); %change in at every time step

% if gexleach <= cexgremain
%     VARIABLES.SOIL.cexgprev = cexgremain-gexleach;
% else
%     VARIABLES.SOIL.cexgprev = 0;
% end
% 
% if fexleach <= cexfremain
%     VARIABLES.SOIL.cexfprev = cexfremain-fexleach;
% else
%     VARIABLES.SOIL.cexfprev = 0;
% end


end