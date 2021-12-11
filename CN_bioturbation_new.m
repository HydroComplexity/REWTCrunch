function [VARIABLES] = CN_bioturbation_new (PARAMS, VARIABLES, CONSTANTS, FORCING, VERTSTRUC, SWITCHES, Cl, Cb, Ch, day)

dz=VERTSTRUC.dzsmm;
CNl=VARIABLES.CNl;
CNl=CNl;
CNh=VARIABLES.Ch./VARIABLES.Nh;
nl_soil=PARAMS.nl_soil;
TBMla=FORCING.TBMla;
TBMla=TBMla(day); % gC/m2

CNveg=PARAMS.CN.CNveg;
CNveg_root=PARAMS.CN.CNveg_root;
% Quijano's code has three crops in one grid. Therefore,
CNveg=mean(CNveg);
CNveg_root=mean(CNveg_root);

%-------------------------------------------------------------------------%
% The same as the main code
%-------------------------------------------------------------------------%


% Total carbon should be inside, not only carbon in litter pool!
Cl_m2=Cl.*(dz./1000);
Ch_m2=Ch.*(dz./1000);
Cb_m2=Cb.*(dz./1000);


% 2 pools
if SWITCHES.CN_type == 0
    Cl_fraction=Cl_m2./sum(Cl_m2+Cb_m2,2);
    Ch_fraction=nan(nl_soil,1);
    Cb_fraction=Cb_m2./sum(Cl_m2+Cb_m2,2);

    C=Cl+Cb;
else % 3 pools
    Cl_fraction=Cl_m2./sum(Cl_m2+Ch_m2+Cb_m2,2);
    Ch_fraction=Ch_m2./sum(Cl_m2+Ch_m2+Cb_m2,2);
    Cb_fraction=Cb_m2./sum(Cl_m2+Ch_m2+Cb_m2,2);

    C=Cl+Ch+Cb;
end


Cflux=TBMla;
fHorO= 0;
% From Cousins et al., 1999
%Dtop=((5*10^-6))*24*365*100*100; % [m2/h] to [cm2/year]
%Dtop=(((5*10^-6)+(1*10^-7))/4)*24*365*100*100; % [m2/h] to [cm2/year]
%Dtop=((10^-7)*5)*24*365*100*100; % [m2/h] to [cm2/year]
%Dtop=((10^-6)/50)*24*365*100*100; % [m2/h] to [cm2/year]
Dtop=((10^-7)*1)*24*365*100*100;%((5*10^-6))*24*365*100*100; % 4 [m2/h] to [cm2/year]
nspecies=1;

% Cumulative depths
zh = cumsum(dz);
% compute nodes for these new layers
zn = (zh + [0 ; zh(1:length(zh)-1)])/2;  

%**************************************************************************
%               BIOTURBATION DUE TO DIFFUSION                             %
%**************************************************************************
BC = zeros(nl_soil,1);
%BC(1) = Cflux*(1-fHorO); % Bioturbation rate 1st layer [grC/m2/dtime]
%BC(2) = Cflux*fHorO;     % Bioturbation rate 2nd layer [grC/m2/dtime]
BC(1) = 0; % Bioturbation rate 1st layer [grC/m2/dtime]
BC(2) = 0;     % Bioturbation rate 2nd layer [grC/m2/dtime]
BC(nl_soil) = 0;

% Bioturbation at the top surface
% Compute Dbio at each interface layer exponential decrease function 
% from Cousins et al., 1999 
D = Dtop*exp(-0.1*(zh/10)); % [cm2/year]  
D = D/100/100/365;          % from [cm2/year] to [m2/day] i.e., [m2/dtime]

% Compute distance between nodes
deltaz = zn(2:nl_soil) - zn(1:nl_soil-1);

%Compute the matrices
[A, B]  = CN_biomatrices  (dz, deltaz, D, 1, BC);
 
% Compute new states of C
%Clnew = A^(-1)*(Cl(1:nl_soil)+B);
%Clnew = Clnew(:);
Cnew = A^(-1)*(C(1:nl_soil)+B);
Cnew = Cnew(:);

% Compute fluxes in and fluxes out in each layer
%[Cin_m2, Cout_m2, difbio_m2, Bioflux] = CN_biofluxes_new (Cnew, dz, deltaz, D, BC);
[Bioflux] = CN_biofluxes_new (Cnew, dz, deltaz, D, BC);

% DK:separte Cl, Ch, and Cb
% Bioflux: 1st colume: up -> down
% Biofulx: 2nd colume: down -> up
if SWITCHES.CN_type == 0
    Cl_Bioflux=Bioflux.*[[0 ; Cl_fraction(1:nl_soil-1)] [Cl_fraction(2:nl_soil); Cl_fraction(nl_soil)]];
    Ch_Bioflux=nan(nl_soil,2);
    Cb_Bioflux=Bioflux.*[[0 ; Cb_fraction(1:nl_soil-1)] [Cb_fraction(2:nl_soil); Cb_fraction(nl_soil)]];
    
    % Put the boundary condition
    Cl_Bioflux(1,1)=Cl_Bioflux(1,1)+Cflux;
else % 3 pools
    Cl_Bioflux=Bioflux.*[[0 ; Cl_fraction(1:nl_soil-1)] [Cl_fraction(1:nl_soil);]];
    Ch_Bioflux=Bioflux.*[[0 ; Ch_fraction(1:nl_soil-1)] [Ch_fraction(1:nl_soil);]];
    Cb_Bioflux=Bioflux.*[[0 ; Cb_fraction(1:nl_soil-1)] [Cb_fraction(1:nl_soil);]];
    
    % Put the boundary condition
    Cl_Bioflux(1,1)=Cl_Bioflux(1,1)+Cflux;
end

clear Bioflux
%%-----------------------------------------------------------------------%%
% Cl case
%difbio_m2 = sum(Bioflux,2);
%Cin_m2_m(Bioflux>0) = Bioflux(Bioflux>0);
%Cout_m2_m(Bioflux<0) = -Bioflux(Bioflux<0);
%Cin_m2 = sum(Cin_m2_m,2);
%Cout_m2 = sum(Cout_m2_m,2);
Cin_m2 = zeros(nl_soil,1);
Cout_m2 = zeros(nl_soil,1);
Cin_m2_m=zeros(nl_soil,2);
Cout_m2_m=zeros(nl_soil,2);

difbio_m2 = sum(Cl_Bioflux,2);
Cin_m2_m(Cl_Bioflux>0) = Cl_Bioflux(Cl_Bioflux>0);
Cout_m2_m(Cl_Bioflux<0) = -Cl_Bioflux(Cl_Bioflux<0);
Cin_m2 = sum(Cin_m2_m,2);
Cout_m2 = sum(Cout_m2_m,2);

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
indmatriz(Cl_Bioflux>0) = Cl_Bioflux(Cl_Bioflux>0);
wgtCinput = (indmatriz./repmat(sum(indmatriz,2),1,2));
indnan = isnan(sum(wgtCinput,2));
% CN 
% matrixinCN = [[CNveg_mean ; CNl(1:nl_soil)] [CNl(1:nl_soil) ; CNl(nl_soil)]];
% matrixinCN = [[CNveg ; CNl(1:nl_soil)] [CNl(1:nl_soil) ; CNl(nl_soil)]];
%matrixinCN = [[CNveg ; CNl(1:nl_soil-1)] [CNl(1:nl_soil-1); CNl(nl_soil)]];
matrixinCN = [[CNveg ; CNl(1:nl_soil-1)] [CNl(2:nl_soil); CNl(nl_soil)]];
%matrixinCN = [[CNveg ; CNl(1:nl_soil-1)] [CNl(1:nl_soil);]];

CN_bio_in = (matrixinCN(:,1).*matrixinCN(:,2))./((matrixinCN(:,2).*wgtCinput(:,1))+((matrixinCN(:,1).*wgtCinput(:,2))));
%CN_bio_in = sum(matrixinCN.*wgtCinput,2);
CN_bio_in(indnan) = CNl(indnan);

% CN
if length(CNl) > 12
    CNl = CNl(2:end);
end
CN_bio_out = CNl;


VARIABLES.CN_bio_in = CN_bio_in;
VARIABLES.CN_bio_out = CN_bio_out;


%%-----------------------------------------------------------------------%%
% Cb case
Cbin_m2 = zeros(nl_soil,1);
Cbout_m2 = zeros(nl_soil,1);
Cbin_m2_m=zeros(nl_soil,2);
Cbout_m2_m=zeros(nl_soil,2);

Cb_difbio_m2 = sum(Cb_Bioflux,2);
Cbin_m2_m(Cb_Bioflux>0) = Cb_Bioflux(Cb_Bioflux>0);
Cbout_m2_m(Cb_Bioflux<0) = -Cb_Bioflux(Cb_Bioflux<0);
Cbin_m2 = sum(Cbin_m2_m,2);
Cbout_m2 = sum(Cbout_m2_m,2);

% Change units to [1/m3]
Cb_difbio_m3 = Cb_difbio_m2./(dz./1000);
Cbin_m3 = Cbin_m2./(dz./1000);
Cbout_m3 = Cbout_m2./(dz./1000);

VARIABLES.Cb_bio_change = Cb_difbio_m3;                % [gr/m3/d]
VARIABLES.Cb_bio_in = Cbin_m3;                        % [gr/m3/d]                    
VARIABLES.Cb_bio_out = Cbout_m3;                      % [gr/m3/d]


if SWITCHES.CN_type == 1
    %%-------------------------------------------------------------------%%
    % Ch case
    Chin_m2 = zeros(nl_soil,1);
    Chout_m2 = zeros(nl_soil,1);
    Chin_m2_m=zeros(nl_soil,2);
    Chout_m2_m=zeros(nl_soil,2);

    Ch_difbio_m2 = sum(Ch_Bioflux,2);
    Chin_m2_m(Ch_Bioflux>0) = Ch_Bioflux(Ch_Bioflux>0);
    Chout_m2_m(Ch_Bioflux<0) = -Ch_Bioflux(Ch_Bioflux<0);
    Chin_m2 = sum(Chin_m2_m,2);
    Chout_m2 = sum(Chout_m2_m,2);
    
    % Change units to [1/m3]
    Ch_difbio_m3 = Ch_difbio_m2./(dz./1000);
    Chin_m3 = Chin_m2./(dz./1000);
    Chout_m3 = Chout_m2./(dz./1000);
    
    % DK: take out due to "no litter pool" and "daily time step"
    % % Expand the matrixes to include litter and change units from [gr/m3/dtime] to [gr/m3/d]
    % difbio_m3 = [ 0 ; difbio_m3].*86400/dtime;
    % Cin_m3 = [ 0 ; Cin_m3].*86400/dtime;
    % Cout_m3 = [ 0 ; Cout_m3].*86400/dtime;
    
    VARIABLES.Ch_bio_change = Ch_difbio_m3;                % [gr/m3/d]
    VARIABLES.Ch_bio_in = Chin_m3;                        % [gr/m3/d]
    VARIABLES.Ch_bio_out = Chout_m3;                      % [gr/m3/d]
    
    % Create weighted fraction of the positive values coming into each layer
    indmatriz = zeros(length(Chin_m3),2);
    indmatriz(Ch_Bioflux>0) = Ch_Bioflux(Ch_Bioflux>0);
    wgtCinput = (indmatriz./repmat(sum(indmatriz,2),1,2));
    indnan = isnan(sum(wgtCinput,2));
    % CN
    % matrixinCN = [[CNveg_mean ; CNl(1:nl_soil)] [CNl(1:nl_soil) ; CNl(nl_soil)]];
    % matrixinCN = [[CNveg ; CNl(1:nl_soil)] [CNl(1:nl_soil) ; CNl(nl_soil)]];
    %matrixinCN = [[CNveg ; CNl(1:nl_soil-1)] [CNl(1:nl_soil-1); CNl(nl_soil)]];
    matrixinCN = [[CNh(1) ; CNh(1:nl_soil-1)] [CNh(2:nl_soil); CNh(nl_soil)]];
    %matrixinCN = [[CNh(1) ; CNh(1:nl_soil-1)] [CNh(1:nl_soil);]];
    
    Ch_CN_bio_in = (matrixinCN(:,1).*matrixinCN(:,2))./((matrixinCN(:,2).*wgtCinput(:,1))+((matrixinCN(:,1).*wgtCinput(:,2))));
    %CN_bio_in = sum(matrixinCN.*wgtCinput,2);
    Ch_CN_bio_in(indnan) = CNh(indnan);
    
    % CN
    Ch_CN_bio_out = CNh;
    
    
    VARIABLES.Ch_CN_bio_in = Ch_CN_bio_in;
    VARIABLES.Ch_CN_bio_out = Ch_CN_bio_out;
    
end

%**************************************************************************
% CHECK MASS BALANCE BY SOLUTION OF BIOTURBATION

% CARBON
%cinput = Cflux.*86400/dtime;     % [grN/m2/d]
cinput = Cflux;     % [grN/m2/d]
%Cchange = (Cin_m3-Cout_m3).*[litterthickness ; dz]; %[grN/m2/d]
% DK: Correct the error!
%Cchange = (Cin_m3-Cout_m3).*dz; %[grN/m2/d]
Cchange = (Cin_m3-Cout_m3).*(dz./1000); %[grN/m2/d]
if SWITCHES.CN_type == 1
    Cchange = ((Cin_m3-Cout_m3)+(Chin_m3-Chout_m3)+(Cbin_m3-Cbout_m3)).*(dz./1000); %[grN/m2/d]
elseif SWITCHES.CN_type == 0
    Cchange = ((Cin_m3-Cout_m3)+(Cbin_m3-Cbout_m3)).*(dz./1000); %[grN/m2/d]
end
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
%Nchange = (Cin_m3./CN_bio_in - Cout_m3./CN_bio_out).*(dz./1000); %[grN/m2/dtime]
if SWITCHES.CN_type == 1
    Nchange = ((Cin_m3./CN_bio_in - Cout_m3./CN_bio_out)+(Chin_m3./Ch_CN_bio_in - Chout_m3./Ch_CN_bio_out)...
        +(Cbin_m3./PARAMS.CN.CNb - Cbout_m3./PARAMS.CN.CNb)).*(dz./1000); %[grN/m2/dtime]
elseif SWITCHES.CN_type == 0
    Nchange = ((Cin_m3./CN_bio_in - Cout_m3./CN_bio_out)...
        +(Cbin_m3./PARAMS.CN.CNb - Cbout_m3./PARAMS.CN.CNb)).*(dz./1000); %[grN/m2/dtime]
end
%VARIABLES.bioNerror = (ninput-sum(Nchange))*86400/dtime;  %[grN/m2/d] 
VARIABLES.bioNerror = (ninput-sum(Nchange));  %[grN/m2/d] 

% SAVE Bioturbation and fragmentation rateS from litter to Horizon.
VARIABLES.Cbiorate = cinput;    
VARIABLES.Nbiorate = ninput;   