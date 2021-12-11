Ka=(1.416+0.01357.*Ts).*(10^(-5));
Kw=(10.^(0.08946+0.03605.*Ts)).*10^(-15);

% NH4m and NH3m are in [mol/L]
syms a;
%for i=1:1:nl_soil
%    b=solve(log10(Ka(i))-log10(Kw(i)) == log10(a)+pH(i), a);
%    NH4m_over_NH3m(i)=b;
%end
%NH4m_over_NH3m=double(NH4m_over_NH3m)';
NH4m_over_NH3m=10.^(log10(Ka)-log10(Kw)-pH);

% The unit of NH4m and NH3m is mol/L
% http://www.convertunits.com/molarmass/NH4
% Molar mass of NH4 = 18.03846 g/mol 
% http://www.convertunits.com/molarmass/NH3
% Molar mass of NH3 = 17.03052 g/mol
NH4_mol_g=1/18.03846;
NH3_mol_g=1/17.03052;

% Find the fraction of NH4m
%for i=1:1:nl_soil
%    b=solve(((NH4_mol_g*(a))/(NH3_mol_g*(1-a))) == NH4m_over_NH3m(i),a);
%    fraction_NH4m(i)=b;
%end
%fraction_NH4m=double(fraction_NH4m)';
f= 1./(NH4m_over_NH3m.*NH3_mol_g./NH4_mol_g+1);
fraction_NH4m=1-f;
fraction_NH4m=fraction_NH4m(:);
% Ensure NH3 amount
for i=1:nl_soil
    if fraction_NH4m(i) >= 0.9999
        fraction_NH4m(i)=0.9999;        
    end
end
    
% NH3 concentration should be some...
for i=1:nl_soil
    if fraction_NH4m(i) >= 0.9999
        fraction_NH4m(i) = 0.9999;
    end
end