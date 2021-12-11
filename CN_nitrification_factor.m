function [F_N] = CN_nitrification_factor(fSn, fTn,fflav)
    
%Nitrif = (1-fflav)'.*(fSn.*fTn*kn.*Cb).* Amm; % gN/m3/d
F_N = (1-fflav)'.*fSn.*fTn;




