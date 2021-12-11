function [] = outtocrunch(primarytocrunch,mintocrunch,soiltemp,F_N,F_D,dz,sm)
% Opens input file and overwrites primary species concentrations, pH, temp conditions

crunchfile = 'ShortCourse1.in';
fid = fopen(crunchfile);
crunchinput = textscan(fid, '%s', 'Delimiter', '\n', 'CollectOutput', true);
fclose(fid);

% Pass soil moisture data
crunchinput{:,:}(18) = strcat('fix_saturation',{' '},num2str(sm));

% Incorporate factors affecting different rates
%Denitrification line (30): denitrification -pathway denitrification  1.0
%-rate 107800 % Susana edit 10/20/2020 %10/21/2020 change to 634 %10/29/20
%changed to 134000
crunchinput{:,:}(31) = strcat('denitrification -pathway denitrification  1.0  -rate',{' '},num2str(F_D*134000));

%Nitrification line (32): Nitrif -pathway Nitrif 1.0 -rate 0.001
%edit: set to zero
crunchinput{:,:}(33) = strcat('Nitrif -pathway Nitrif 1.0 -rate',{' '},num2str(F_N*0.001));

%Update temperature for each layer
crunchinput{:,:}{38} = soiltemp; 

%Overwrite initial condition block (primary species):
crunchinput{:,:}(39:56) = primarytocrunch;
crunchinput{:,:}(58:71) = mintocrunch;

%Overwrite layer depth (cm)
crunchinput{:,:}(85) = strcat('xzones 1',{' '},num2str(dz.*100)); 

newinput = 'ShortCourse1.in';
dlmcell(newinput,crunchinput{:})
fclose('all');

end