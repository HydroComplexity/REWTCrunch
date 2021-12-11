%*************************************************************************
%                   COMPUTE ENVIRONMENTAL FACTORS
%*************************************************************************
    
indT = Ts > 25;
fTd(indT) = 1;              
fTd(~indT) = 1.9.^((Ts(~indT)-25)/10); 
fTd = fTd(:);

% Calculate soil matric potential from soil moisture
for i=1:nl_soil
   if soiltype(i) == 9203 % sand
       res_sm=0.058; % [cm3/cm3]
%       sat_sm=0.37; % [cm3/cm3]
       sat_sm=0.4370001; % [cm3/cm3]
       alpha=0.035;  % [1/cm]
       n=3.19;
       BulkDens=1.69; % [g/cm3]
       FieldCapa=0.12; % [-]
       ClayContent=0.03;
   end
   if soiltype(i) == 8107 % loamy sand
       res_sm=0.074;
%       sat_sm=0.39;
       sat_sm=0.4370001;
       alpha=0.035;
       n=2.39;
       BulkDens=1.63; % [g/cm3]
       FieldCapa=0.14; % [-]
       ClayContent=0.07;
   end
   if soiltype(i) == 6510 % sandy loam
       res_sm=0.067;
%       sat_sm=0.37;
       sat_sm=0.4530001;
       alpha=0.021;
       n=1.61;
       BulkDens=1.51; % [g/cm3]
       FieldCapa=0.23; % [-]
       ClayContent=0.10;
   end
   if soiltype(i) == 4218 % loam
       res_sm=0.083;
       %sat_sm=0.46;
       sat_sm=0.463;
       alpha=0.025;
       n=1.31;
       BulkDens=1.43; % [g/cm3]
       FieldCapa=0.26; % [-]
       ClayContent=0.18;
   end
   
   if soiltype(i) == 2015 % silt loam
       res_sm=0.061;
%       sat_sm=0.43;
       sat_sm= 0.501;
       alpha=0.012;
       n=1.39;
       BulkDens=1.38; % [g/cm3]
       FieldCapa=0.3; % [-]
       ClayContent=0.15;
   end
   if soiltype(i) == 6027 % sandy clay loam
       res_sm=0.086;
%       sat_sm=0.40;
       sat_sm=0.3980001;
       alpha=0.033;
       n=1.49;
       BulkDens=1.40; % [g/cm3]
       FieldCapa=0.33; % [-]
       ClayContent=0.27;
   end
   if soiltype(i) == 3234 % clay loam
       res_sm=0.129;
%       sat_sm=0.47;
       sat_sm=0.464;
       alpha=0.030;
       n=1.37;
       BulkDens=1.30; % [g/cm3]
       FieldCapa=0.335; % [-]
       ClayContent=0.34;
   end
   if soiltype(i) == 933 % silty clay loam
       res_sm=0.098;
%       sat_sm=0.55;
       sat_sm=0.471;
       alpha=0.027;
       n=1.41;
       BulkDens=1.27; % [g/cm3]
       FieldCapa=0.34; % [-]
       ClayContent=0.33;
   end
   if soiltype(i) == 1045 % silty clay
       res_sm=0.163;
%       sat_sm=0.47;
       sat_sm=0.4790001;
       alpha=0.023;
       n=1.39;
       BulkDens=1.21; % [g/cm3]
       FieldCapa=0.36; % [-]
       ClayContent=0.45;
   end
   if soiltype(i) == 2060 % clay
       res_sm=0.102;
%       sat_sm=0.51;
       sat_sm=0.4750001;
       alpha=0.021;
       n=1.20;
       BulkDens=1.25; % [g/cm3]
       FieldCapa=0.36; % [-]
       ClayContent=0.60;
   end
   
   % For matric_potential
   if real_sm(i) < res_sm
       real_sm(i)=res_sm+0.000000000000001;
   end
   matric_potential(i)=(1/alpha)*((1/(((real_sm(i)-res_sm)/(sat_sm-res_sm))^(1/(1-1/n))))-1)^(1/n); % [cm]
   if isreal(matric_potential(i)) ==0
      stop=1; 
   end
   
   % Water-Filled Pore Space
   WFPS(i)=(sm(i)./sat_sm); 
   indOver1=WFPS>1;
   WFPS(indOver1)=1;
   
   % Clay content
   Clay(i)=ClayContent;
   
   %matric_potential(i)=(1/alpha)*((1/(((real_sm(i)-res_sm)/(sat_sm-res_sm))^(1/(1-1/n))))-1)^(1/n); % [cm]
   BulkDensity(i)=BulkDens.*(100*100*100); % [g/cm3] -> [g/m3]
   BulkDensity=BulkDensity(:);
   
   % porosity
   porosity(i)=sat_sm;
   
   FieldCapacity(i)=FieldCapa;
   FieldCapacity=FieldCapacity(:);
   
   %SoilAir(i)=porosity(i)-sm(i);
   %if SoilAir(i)<0
   %    SoilAir(i)=0;
   %end
   %Dfc(i)=porosity(i)^2*((SoilAir(i))/porosity(i))^(2+3/(13.6*ClayContent+3.5));
   %Dfc=Dfc(:);

end

PARAMS.porosity=porosity';
WFPS=WFPS';

% http://www.convertunits.com/from/cm+H2O/to/megapascal
smp=matric_potential.*0.0000980665; % [cm] to [MPa]
fSd = log(7.58./abs(smp))/log(7.58/0.01);
indWP = fSd < 0;
fSd(indWP) = 0;
indWP = fSd > 1;
fSd(indWP) = 1;
fSd = fSd(:);

% Porporato 2003
% for i=1:nl_soil
%    if sm(i) <= FieldCapacity(i)
%        fSd(i)=sm(i)/FieldCapacity(i);
%    else
%        fSd(i)=FieldCapacity(i)/sm(i);       
%    end
% end
%fSd(17:end)=0;

% Wrong!
%smp = smp * 9.8066e-06; % First convert to MPa

% DK: fTd, Paul, K.(2001)
%  / 1/[1+ 4. exp(-6.RWC)] / water content relative to 
% the upper and lower limit of water observed
%in the field (i.e. RWC=(W-LL)/(UL-LL)).
% for iiii=1:nl_soil
%     fSd(iiii) = 1/(1+4*exp(-6*((real_sm(iiii)-0)/(porosity(iiii)-0)))); 
% end
% fSd = fSd(:);


% for Nitrification
fTn = fTd;
fSn = fSd;

% Porporato 2003
% for i=1:nl_soil
%    if sm(i) <= FieldCapacity(i)
%        fSn(i)=sm(i)/FieldCapacity(i);
%    else
%        fSn(i)=(1-sm(i))/(1-FieldCapacity(i));       
%    end
% end
porosity=porosity(:);
if SWITCHES.CN.PoporatofSdfNd_on == 1
    indS = real_sm(1:end-1) > FieldCapacity;
    fSd(indS) = (FieldCapacity(indS)./porosity(indS))./(real_sm(indS)./porosity(indS));
    fSd(~indS) = (real_sm(~indS)./porosity(~indS))./(FieldCapacity(~indS)./porosity(~indS));
    fSd=fSd(:);
    
    indSn = real_sm(1:end-1) > FieldCapacity;
    fSn(indSn) = (1-(real_sm(indSn)./porosity(indSn)))./(1-(FieldCapacity(indSn)./porosity(indSn)));
    fSn(~indSn) = (real_sm(~indSn)./porosity(~indSn))./(FieldCapacity(~indSn)./porosity(~indSn));
    fSn=fSn(:);   
    fSn(fSn<0.000001)=0;
end
