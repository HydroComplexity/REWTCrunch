function [ fflav ] = CN_BNIfactor( i1,i2,cexf)

Asoil = 1; %sq unit
 
% cexf = cexf/270.24/1000; % mmol/ g dry soil/day (MW_genistein = 270.24 g/mol (Pubchem))

 for i = i1:i2
% Need g of soil into which flavonoids are released. 
%     exflav = cexf>0;
%     gsoil = sum(Asoil.*dz(exflav).*BulkDensity(exflav));
    ATU(i) = cexf(i)/(10^-6)*10; %conversion from sorghum study with sakuranetin
%     fflav = 1/100*(23.01*log(ATU)+0.4735); %fraction, not percent
%     fflav = 1/100*(4.7*ATU - 10.6);
    if ATU(i) > 0.0001
        fflav(i) = 1/100*(15.948*ATU(i)^0.5174);
%         fflav = 1/100*(-0.2101*ATU^2 + 6.8512*ATU + 7.727);
    else
        fflav(i) = 0;
    end
    if fflav(i) < 0
        fflav(i) = 0;
    end
    if fflav(i) > 1
        fflav(i) = 1;
    end
%     if isinf(fflav(i))
%         fflav(i) = 0;
%     elseif fflav(i) <= 1
%         fflav(i) = 0;        
%     end
 end


end

