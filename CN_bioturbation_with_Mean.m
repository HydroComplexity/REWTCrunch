function [VARIABLES] = CN_bioturbation_with_Mean (PARAMS, VARIABLES...
    , CONSTANTS, FORCING, VERTSTRUC, Cl, TBMla, CNvegtation, CNvegtation_root, Clpropotion)

dz=VERTSTRUC.dzsmm;
CNl=VARIABLES.CNl;
CNl = CNl(2:end);
nl_soil=PARAMS.nl_soil;
% Put average value
%TBMla=TBMla;
%TBMla=TBMla(day);
TBMla=TBMla;

% Put average value
%CNveg=PARAMS.CN.CNveg;
CNveg=CNvegtation;
CNveg_root=CNvegtation_root;

% Quijano's code has three crops in one grid. Therefore,
CNveg=mean(CNveg);
CNveg_root=mean(CNveg_root);

Cflux=TBMla;
fHorO= 0;
% From Cousins et al., 1999
%Dtop=((5*10^-6))*24*365*100*100; % [m2/h] to [cm2/year]
%Dtop=(((5*10^-6)+(1*10^-7))/2)*24*365*100*100; % [m2/h] to [cm2/year]
%Dtop=((10^-7)*5)*24*365*100*100; % [m2/h] to [cm2/year]
%Dtop=((10^-6)/50)*24*365*100*100; % [m2/h] to [cm2/year] = 1.75
%Dtop=((10^-6)/50)*24*365*100*100; % [m2/h] to [cm2/year]
Dtop=((10^-7)*1)*24*365*100*100; %4
nspecies=1;

% Cumulative depths
zh = cumsum(dz);
% compute nodes for these new layers
zn = (zh + [0 ; zh(1:length(zh)-1)])/2;  

%**************************************************************************
%               BIOTURBATION DUE TO DIFFUSION                             %
%**************************************************************************
BC = zeros(nl_soil,1);
BC(1) = Cflux*(1-fHorO); % Bioturbation rate 1st layer [grC/m2/dtime]
BC(2) = Cflux*fHorO;     % Bioturbation rate 2nd layer [grC/m2/dtime]
BC(nl_soil) = 0;

% Bioturbation at the top surface
% Compute Dbio at each interface layer exponential decrease function 
% from Cousins et al., 1999 
D = Dtop*exp(-0.1*(zh/10)); % [cm2/year]  
D = D/100/100/365;          % from [cm2/year] to [m2/day] i.e., [m2/dtime]

% Compute distance between nodes
deltaz = zn(2:nl_soil) - zn(1:nl_soil-1);

%Compute the matrices
[A, B]  = CN_biomatrices_with_Mean (dz, deltaz, D, 1, BC);
 
% Compute new states of C
Clnew = A^(-1)*(Cl(1:nl_soil)+B);
Clnew = Clnew(:);

% Compute fluxes in and fluxes out in each layer
[Cin_m2, Cout_m2, difbio_m2, Bioflux] = CN_biofluxes_with_Mean (Clnew, dz, deltaz, D, BC);

% Change units to [1/m3]
difbio_m3 = difbio_m2./(dz./1000);
Cin_m3 = Cin_m2./(dz./1000);
Cout_m3 = Cout_m2./(dz./1000);

% DK: take out due to "no litter pool" and "daily time step"
% % Expand the matrixes to include litter and change units from [gr/m3/dtime] to [gr/m3/d]
% difbio_m3 = [ 0 ; difbio_m3].*86400/dtime;
% Cin_m3 = [ 0 ; Cin_m3].*86400/dtime;
% Cout_m3 = [ 0 ; Cout_m3].*86400/dtime;

VARIABLES.Cl_bio_change = difbio_m3;                % [gr/m3/d]
VARIABLES.Cl_bio_in = Cin_m3;                        % [gr/m3/d]                    
VARIABLES.Cl_bio_out = Cout_m3;                      % [gr/m3/d]

% Create weighted fraction of the positive values coming into each layer 
indmatriz = zeros(length(Cin_m3),2);
indmatriz(Bioflux>0) = Bioflux(Bioflux>0);
wgtCinput = (indmatriz./repmat(sum(indmatriz,2),1,2));
indnan = isnan(sum(wgtCinput,2));
% CN 
% matrixinCN = [[CNveg_mean ; CNl(1:nl_soil)] [CNl(1:nl_soil) ; CNl(nl_soil)]];
% matrixinCN = [[CNveg ; CNl(1:nl_soil)] [CNl(1:nl_soil) ; CNl(nl_soil)]];
%matrixinCN = [[CNveg ; CNl(1:nl_soil-1)] [CNl(1:nl_soil-1); CNl(nl_soil)]];
matrixinCN = [[CNveg ; CNl(1:nl_soil-1)] [CNl(2:nl_soil); CNl(nl_soil)]];

CN_bio_in = (matrixinCN(:,1).*matrixinCN(:,2))./((matrixinCN(:,2).*wgtCinput(:,1))+((matrixinCN(:,1).*wgtCinput(:,2))));

%CN_bio_in = sum(matrixinCN.*wgtCinput,2);
CN_bio_in(indnan) = CNl(indnan);
% CN
CN_bio_out = CNl;
VARIABLES.CN_bio_in = CN_bio_in;
VARIABLES.CN_bio_out = CN_bio_out;

%**************************************************************************
% CHECK MASS BALANCE BY SOLUTION OF BIOTURBATION

% CARBON
%cinput = Cflux.*86400/dtime;     % [grN/m2/d]
cinput = Cflux;     % [grN/m2/d]
%Cchange = (Cin_m3-Cout_m3).*[litterthickness ; dz]; %[grN/m2/d]
% DK: Correct the error!
%Cchange = (Cin_m3-Cout_m3).*dz; %[grN/m2/d]
Cchange = (Cin_m3-Cout_m3).*(dz./1000); %[grN/m2/d]
VARIABLES.bioCerror = (cinput-sum(Cchange));  %[grC/m2/d] 

% NITROGEN
%ninput = Cflux/(CNl(1)).*86400/dtime;    % [grN/m2/dtime]
%ninput = Cflux/(CNl(1));    % [grN/m2/dtime]
ninput = Cflux/(CNveg);    % [grN/m2/dtime]
%Nin = (Cin_m3.*[litterthickness ; dz])./CN_bio_in;
% DK: error corection
%Nin = (Cin_m3.*dz)./CN_bio_in;
Nin = (Cin_m3.*(dz./1000))./CN_bio_in;
%Nout = (Cout_m3.*[litterthickness ; dz])./CN_bio_out;
% DK: error corection
%Nout = (Cout_m3.*dz)./CN_bio_out;
%Nchange = (Cin_m3./CN_bio_in - Cout_m3./CN_bio_out).*dz; %[grN/m2/dtime]
Nout = (Cout_m3.*(dz./1000))./CN_bio_out;
Nchange = (Cin_m3./CN_bio_in - Cout_m3./CN_bio_out).*(dz./1000); %[grN/m2/dtime]
%VARIABLES.bioNerror = (ninput-sum(Nchange))*86400/dtime;  %[grN/m2/d] 
VARIABLES.bioNerror = (ninput-sum(Nchange));  %[grN/m2/d] 

% SAVE Bioturbation and fragmentation rateS from litter to Horizon.
VARIABLES.Cbiorate = cinput;    
VARIABLES.Nbiorate = ninput;   