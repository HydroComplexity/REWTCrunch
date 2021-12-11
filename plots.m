%% Plotting script for REWTCrunch
% S. Roque-Malo

close all
clear all
clc

% Open REWTCrunch and REWT results files
crex = load('your_REWTCrunch_with_exudation_results');
crnoex = load('your_REWTCrunch_withOUT_exudation_results');
REWT = load('your_REWT_withexudation_results');
REWTnoex = load('your_REWT_withOUT_exudation_results');


%Other necessary variables
z95 = 0.9;
botrow = 12;
maxrf = 0.2;
spec = 'corn';
deltick = [21:50:171];
dat = 'del';

% Plotting stuff
 XX = 1:length(crex.Nit_print(1,:));
 XXr = 1:length(crex.Nit_print(1,:));
 Y = crex.Nit; 
 Z = NaN(size(Y));
 ZZ = cumsum(crex.dz(1:botrow));
 Z = repmat(ZZ,1,max(XX));
 X = repmat(XX,botrow,1); 
 Xr = repmat(XXr,botrow,1);
 depths = crex.depths(1:botrow);
%  

% Initial values (These must match the values in main code
SiO2_o = 4.18006641001185E-04;
Al3p_o = 5.52480585268130E-08;
Mg2p_o = 1.96229806315877E-04;
Ca2p_o = 1.38292166758517E-03;
Kp_o = 2.22547613802969E-06;
Fe2p_o = 7.55E-08;
Fe3p_o = 6.39E-22;
TiO4H4_o = 2.76081894084137E-09;
Nap_o = 7.12193976883976E-04;
NH4_o = 0.52;
Nit_o = 2;

%% Plots

% Ammonia, Nitrate, Ammonium concentrations (Figure 4)
sm=crex.vwc_layer_day_point(1:end-1,:)+crex.ice_layer_day_point(1:end-1,:);
figure()
h(1) = subplot(3,12,1:5)    ;
    crex.Amm_print = crex.Amm_print./18.03846.*0.001./sm;
    [datlargestcr,datsmallestcr,newdat] = notails(crex.Amm_print,100);
    pcolor(Xr,ZZ,log10(crex.Amm_print));
    alcb = colorbar;
    alcb.Ticks = [ -25 -15 -5 ];
    alcb.TickLabels = {'10^{-25}' '10^{-15}' '10^{-5}' };
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('REWTCrunch - NH_4^+')

 h(2) = subplot(3,12,7:11)    ;
    REWT.Amm_print = REWT.Amm_print./18.03846.*0.001./sm;
    [datlargest,datsmallest,newdat] = notails(REWT.Amm_print,100);
    pcolor(Xr,ZZ,log10(REWT.Amm_print));
    alcb = colorbar;
    alcb.Ticks = [-8 -6 -4 -2  ];
    alcb.TickLabels = {'10^{-8}' '10^{-6}' '10^{-4}' '10^{-2}' };
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('REWT - NH_4^+')
    
 h(3) = subplot(3,12,13:17)    ;
    crex.Nit_print = crex.Nit_print./62.0049.*0.001./sm;
    [datlargestcr,datsmallestcr,newdat] = notails(crex.Nit_print,10);
    pcolor(Xr,ZZ,log10(crex.Nit_print));
    alcb = colorbar;
    alcb.Ticks = [-12 -8 -4   ];
    alcb.TickLabels = {'10^{-12}' '10^{-8}' '10^{-4}' };
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('REWTCrunch - NO_3^-')
    
h(4) = subplot(3,12,19:23)    ;
    REWT.Nit_print = REWT.Nit_print./62.0049.*0.001./sm;
    [datlargest,datsmallest,newdat] = notails(REWT.Nit_print,100);
    pcolor(Xr,ZZ,log10(REWT.Nit_print));
    alcb = colorbar;
    alcb.Ticks = [-7 -5 -3 ];
    alcb.TickLabels = {'10^{-7}' '10^{-5}' '10^{-3}' };
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('REWT - NO_3^-')
    
 h(5) = subplot(3,12,25:29)    ;
    crex.NH3_print = crex.NH3_print./17.03052.*0.001./sm;
    [datlargestcr,datsmallestcr,~] = notails(crex.NH3_print,100);
    pcolor(Xr,ZZ,log10(crex.NH3_print));
    alcb = colorbar;
    alcb.Ticks = [-6 -5 -4  -3 ];
    alcb.TickLabels = {'10^{-6}' '10^{-5}' '10^{-4}' '10^{-3}'  };
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('REWTCrunch - NH_3')
    xlabel('Days of Year')
    
h(6) = subplot(3,12,31:35)    ;
    REWT.NH3_print = REWT.NH3_print./17.03052.*0.001./sm;
    [datlargest,datsmallest,newdat] = notails(REWT.NH3_print,100);
    pcolor(Xr,ZZ,log10(REWT.NH3_print));
    alcb = colorbar;
    alcb.Ticks = [-20 -15 -10 ];
    alcb.TickLabels = {'10^{-20}' '10^{-15}' '10^{-10}' };
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('REWT - NH_3')
    xlabel('Days of Year')
set(findall(gcf,'-property','FontSize'),'FontSize',12)

% % Flavonoid concentration (mol/kg water) (Figure 5) 

figure()
h(1) = subplot(1,2,1)    ;
    crex.cexf_print = crex.cexf_print./270.24.*0.001./sm;
    [datlargest,datsmallest,newdat] = notails(crex.cexf_print,95);
    pcolor(Xr,ZZ,crex.cexf_print);
    caxis([datsmallest datlargest])
    alcb = colorbar;
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlabel('DoY')
    xlim([1 365])
    xlabel('Days of Year')
    title('REWTCrunch')
set(gca, 'FontSize', 12);
 h(2) = subplot(1,2,2)    ;
    REWT.cexf_print = REWT.cexf_print./270.24.*0.001./sm;
    [datlargest,datsmallest,newdat] = notails(REWT.cexf_print,95);
    pcolor(Xr,ZZ,REWT.cexf_print);
    alcb = colorbar;
    caxis([datsmallest datlargest])
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlabel('DoY')
    xlim([1 365])
    xlabel('Days of Year')
    title('REWT')
   set(gca, 'FontSize', 12);

% Microbial biomass (gC/m3) (Figure 6)
figure()
h1 = subplot(1,2,1);
    [datlargestcr,datsmallestcr,newdat] = notails(crex.Cb_print,100);
    pcolor(Xr,ZZ,crex.Cb_print);
    colorbar
    caxis([datsmallestcr datlargestcr])
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('REWTCrunch - Microbial Biomass')
    xlabel('Day of Year')
h(2)=subplot(1,2,2)
    pchcb = (crex.Cb_print-REWT.Cb_print)./REWT.Cb_print.*100;
    customCMap1 = customcolormap(pchcb,0,max(max(pchcb)),min(min(pchcb)));    
    pcolor(Xr,ZZ,pchcb);
    colormap(h(2),customCMap1)
    colorbar
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(%)')
    xlim([1 365])
    title({'Percent Change';' '})
    xlabel('Day of Year')
set(findall(gcf,'-property','FontSize'),'FontSize',12)

% SiO2 figure (Figure 7)
figure()
h(1) = subplot(2,11,1:5)    ;
    [datlargest,datsmallest,newdat] = notails(crex.SiO2_print,100);
    pcolor(Xr,ZZ,log10(crex.SiO2_print));
    alcb = colorbar;
    alcb.Ticks = [-5 -4.5 -4 -3.5 -3];
    alcb.TickLabels = {'10^{-5}' '10^{-4.5}' '10^{-4}' '10^{-3.5}' '10^{-3}'};
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlabel('DoY')
    xlim([1 365])
    title({'SiO_2_(_a_q_)';'(mol/kg water)'})
h(2) = subplot(2,11,7:11)    ;
    pltdif = (crex.SiO2_print-crnoex.SiO2_print)./crnoex.SiO2_print.*100;
    pltdif(isnan(pltdif))=0;
    pltdif(isinf(pltdif))=0;
    pltdif = sign(pltdif).*log10(1+abs(pltdif)/10^-2);
    [datlargest,datsmallest,newdat] = notails(pltdif,100);
    customCMap1 = customcolormap(newdat,0,datlargest,datsmallest);
    pcolor(Xr,ZZ,pltdif);
    colormap(h(2),customCMap1)
    cbh = colorbar;
    cbh.Ticks = [-3 -2 -1 0 1 2 3];
    cbh.TickLabels = num2cell([-10 -1 -0.1 0 0.1 1 10]);
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlabel('DoY')
    xlim([1 365])
    title('Change in SiO_2_(_a_q_) due to Root Exudation')
h(4) = subplot(2,11,12:16)    ;
    discharge = crex.qq_layer_day_point(13,:);
    janmar = 1:90;
    aprjun = 91:181;
    julsep = 182:273;
    octdec = 274:365;
    loglog(discharge(janmar),crex.TLCH_SiO2_m2_print(janmar),'.')
    hold on;
    loglog(discharge(aprjun),crex.TLCH_SiO2_m2_print(aprjun),'.')
    hold on;
    loglog(discharge(julsep),crex.TLCH_SiO2_m2_print(julsep),'.')
    hold on;
    loglog(discharge(octdec),crex.TLCH_SiO2_m2_print(octdec),'.')
    title('SiO_2_(_a_q_) Leaching Concentration vs. Leaching Flux')
    xlabel('Log(Flux)')
    ylabel('Log(Concentration)')
    legend([h(4)],{'Jan-Mar','Apr-Jun','Jul-Sep','Oct-Dec'},'Location','NorthWest');
h(3) =  subplot(2,11,18:22)     ;
    pltdif = (crex.TLCH_SiO2_m2_print-crnoex.TLCH_SiO2_m2_print)./crnoex.TLCH_SiO2_m2_print.*100;
    pltdif(isnan(pltdif))=0;
    pltdif(isinf(pltdif))=0;
    plot(pltdif)
    colorbar
    ylabel('(%)')
    xlim([1 365])
    xlabel('DoY')
    title('Change in SiO_2_(_a_q_) Leaching due to Root Exudation')
set(findall(gcf,'-property','FontSize'),'FontSize',12)

%%
% Ca2+ figure (Figure 8)
figure()
h(1) = subplot(2,11,1:5)    ;
    [datlargest,datsmallest,newdat] = notails(crex.Ca2p_print,100);
    pcolor(Xr,ZZ,log10(crex.Ca2p_print));
    alcb = colorbar;
    alcb.Ticks = [-6 -5 -4 -3];
    alcb.TickLabels = {'10^{-6}' '10^{-5}' '10^{-4}' '10^{-3}'};
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlabel('DoY')
    xlim([1 365])
    title('Ca^2^+ [mol/kg water]')

h(2) = subplot(2,11,7:11)    ;
    pltdif = (crex.Ca2p_print-crnoex.Ca2p_print)./crnoex.Ca2p_print.*100;
    pltdif(isnan(pltdif))=0;
    pltdif(isinf(pltdif))=0;
     pltdif = sign(pltdif).*log10(1+abs(pltdif)/10^-2);
    [datlargest,datsmallest,newdat] = notails(pltdif,100);
    customCMap1 = customcolormap(newdat,0,datlargest,datsmallest);
    pcolor(Xr,ZZ,pltdif);
    colormap(h(2),customCMap1)
    cbh = colorbar;
    cbh.Ticks = [-1 0 1 2 3];
    cbh.TickLabels = num2cell([-0.1 0 0.1 1 10]);
     shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    xlabel('DoY')
    title('Change in Ca^2^+ due to Root Exudation')
    
h(4) = subplot(2,11,12:16)   ;
    discharge = crex.qq_layer_day_point(13,:);
    janmar = 1:90;
    aprjun = 91:181;
    julsep = 182:273;
    octdec = 274:365;
    loglog(discharge(janmar),crex.TLCH_Ca2p_m2_print(janmar),'.')
    hold on;
    loglog(discharge(aprjun),crex.TLCH_Ca2p_m2_print(aprjun),'.')
    hold on;
    loglog(discharge(julsep),crex.TLCH_Ca2p_m2_print(julsep),'.')
    hold on;
    loglog(discharge(octdec),crex.TLCH_Ca2p_m2_print(octdec),'.')
    colorbar
    title('Ca^2^+ Leaching Concentration vs. Leaching Flux')
    xlabel('Log(Flux)')
    ylabel('Log(Concentration)')
    legend([h(4)],{'Jan-Mar','Apr-Jun','Jul-Sep','Oct-Dec'},'Location','NorthWest');
h(3) =    subplot(2,11,18:22) ;
    pltdif = (crex.TLCH_Ca2p_m2_print-crnoex.TLCH_Ca2p_m2_print)./crnoex.TLCH_Ca2p_m2_print.*100;
    pltdif(isnan(pltdif))=0;
    pltdif(isinf(pltdif))=0;
    plot(pltdif)
    colorbar
    ylabel('(%)')
    xlim([1 365])
    xlabel('DoY')
    title('Change in Ca^2^+ Leaching due to Root Exudation')
    
set(findall(gcf,'-property','FontSize'),'FontSize',12)

%%
% Al3+ figure (Figure 9)
figure()
h(1) = subplot(2,11,1:5)    ;
    [datlargest,datsmallest,newdat] = notails(crex.Al3p_print,100);
    pcolor(Xr,ZZ,log10(crex.Al3p_print));
    alcb = colorbar;
    alcb.Ticks = [-12 -10 -8 -6 -4];
    alcb.TickLabels = {'10^{-12}' '10^{-10}' '10^{-8}' '10^{-6}' '10^{-4}'};
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlabel('DoY')
    xlim([1 365])
    title({'Al^3^+';'(mol/kg water)'})
h(2) = subplot(2,11,7:11)    ;
    pltdif = (crex.Al3p_print-crnoex.Al3p_print)./crnoex.Al3p_print.*100;
    pltdif(isnan(pltdif))=0;
    pltdif(isinf(pltdif))=0;    
    pltdif = sign(pltdif).*log10(1+abs(pltdif)/10^-2);
    [datlargest,datsmallest,newdat] = notails(pltdif,100);
    customCMap1 = customcolormap(newdat,0,datlargest,datsmallest);
    pcolor(Xr,ZZ,pltdif);
    colormap(h(2),customCMap1)
    cbh = colorbar;
    cbh.Ticks = [-2 0 2 4 6];
    cbh.TickLabels = num2cell([-1 0 1 100 1000]);
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlabel('DoY')
    xlim([1 365])
    title('Change in Al^3^+ due to Root Exudation')
 h(4) = subplot(2,11,12:16)     ;
    discharge = crex.qq_layer_day_point(13,:);
    janmar = 1:90;
    aprjun = 91:181;
    julsep = 182:273;
    octdec = 274:365;
    loglog(discharge(janmar),crex.TLCH_Al3p_m2_print(janmar),'.')
    hold on;
    loglog(discharge(aprjun),crex.TLCH_Al3p_m2_print(aprjun),'.')
    hold on;
    loglog(discharge(julsep),crex.TLCH_Al3p_m2_print(julsep),'.')
    hold on;
    loglog(discharge(octdec),crex.TLCH_Al3p_m2_print(octdec),'.')
    colorbar
    title('Al^3^+ Leaching Concentration vs. Leaching Flux')
    xlabel('Log(Flux)')
    ylabel('Log(Concentration)')
    legend([h(4)],{'Jan-Mar','Apr-Jun','Jul-Sep','Oct-Dec'},'Location','NorthWest');
h(3) =  subplot(2,11,18:22)   ;
    pltdif = (crex.TLCH_Al3p_m2_print-crnoex.TLCH_Al3p_m2_print)./crnoex.TLCH_Al3p_m2_print.*100;
    pltdif(isnan(pltdif))=0;
    pltdif(isinf(pltdif))=0;
    plot(pltdif)
    colorbar
    ylabel('(%)')
    xlim([1 365])
    title('Change in Al^3^+ Leaching due to Root Exudation')
    
set(findall(gcf,'-property','FontSize'),'FontSize',12)


% Nitrate leaching comparison (Figure 10)
figure()
h1 = subplot(2,1,1);
    RNL = sign(REWT.TLCH_nit_m2_print).*log10(1+abs(REWT.TLCH_nit_m2_print)/10^-4);
    CRNL = sign(crex.TLCH_Nit_m2_print).*log10(1+abs(crex.TLCH_Nit_m2_print)/10^-4);
    plot(RNL,'+')
    hold on;
    plot(CRNL,'.') 
    ylabel('(g/m^2)')
    xlim([1 365])
    ax = gca();
    yticks([-1 0 1 2 3]);
    yticklabels([{'-0.001' '0' '0.001' '0.01' '0.1'}]);
    legend([h1],{'REWT - exudation','REWTCrunch - exudation'},'Location','NorthEast');
    title('Total Nitrate Leached')
h2 = subplot(2,1,2);
    rewtexdiff = (REWT.TLCH_nit_m2_print-REWTnoex.TLCH_nit_m2_print)./(REWTnoex.TLCH_nit_m2_print).*100;
    crrewtnoexdiff = (crnoex.TLCH_Nit_m2_print-REWTnoex.TLCH_nit_m2_print)./REWTnoex.TLCH_nit_m2_print.*100;
    crrewtexdiff = (crex.TLCH_Nit_m2_print-REWTnoex.TLCH_nit_m2_print)./REWTnoex.TLCH_nit_m2_print.*100;
    plot(crrewtnoexdiff);
    hold on;
    plot(crrewtexdiff);
    legend([h2],{'REWTCrunch - no exudation','REWTCrunch - exudation'},'Location','NorthEast');
    xlim([1 365])
    xlabel('Days of Year')
    ylabel('(%)')
    title({'Percent Change Compared to "REWT-no exudation"';' '})
set(findall(gcf,'-property','FontSize'),'FontSize',12)



% Mineral solute concentrations (mol/kg water) (Figure 7, 8, and 9 (a) and Figure A.3)
figure()
h(1) = subplot(5,2,1)    ;
    [datlargest,datsmallest,newdat] = notails(crex.SiO2_print,100);
    pcolor(Xr,ZZ,log10(crex.SiO2_print));
    alcb = colorbar;
    alcb.Ticks = [-5 -4 -3];
    alcb.TickLabels = {'10^{-5}' '10^{-4}' '10^{-3}'};
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('SiO_2_(_a_q_)')

h(2) = subplot(5,2,2) ;
    [datlargest,datsmallest,newdat] = notails(crex.TiO4H4_print,100);
    pcolor(Xr,ZZ,log10(crex.TiO4H4_print));
    alcb = colorbar;
    alcb.Ticks = [ -9.301 -9 -8.698 -8.52];
    alcb.TickLabels = {'0.5' '1' '2' '3'};
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('Ti(OH)_4_(_a_q_)')

h(3) = subplot(5,2,3)  ;
    [datlargest,datsmallest,newdat] = notails(crex.Ca2p_print,100);
    pcolor(Xr,ZZ,log10(crex.Ca2p_print));
    alcb = colorbar;
    alcb.Ticks = [-6 -5 -4 -3 ];
    alcb.TickLabels = {'10^{-6}' '10^{-5}' '10^{-4}' '10^{-3}' };
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('Ca^2^+')

h(4) = subplot(5,2,4);
    [datlargest,datsmallest,newdat] = notails(crex.Mg2p_print,100);
    pcolor(Xr,ZZ,log10(crex.Mg2p_print));
    alcb = colorbar;
    alcb.Ticks = [-8 -6 -4];
    alcb.TickLabels = {'10^{-8}' '10^{-6}' '10^{-4}'};
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('Mg^2^+')

h(5) = subplot(5,2,5)  ;
    [datlargest,datsmallest,newdat] = notails(crex.Kp_print,100);
    pcolor(Xr,ZZ,log10(crex.Kp_print));
    alcb = colorbar;
    alcb.Ticks = [-7 -6 -5 -4 -3];
    alcb.TickLabels = {'10^{-7}' '10^{-6}' '10^{-5}' '10^{-4}' '10^{-3}'};
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('K^+')

h(6) = subplot(5,2,6)  ;
    [datlargest,datsmallest,newdat] = notails(crex.Nap_print,100);
    pcolor(Xr,ZZ,log10(crex.Nap_print));
    alcb = colorbar;
    alcb.Ticks = [-7 -6 -5 -4 -3];
    alcb.TickLabels = {'10^{-7}' '10^{-6}' '10^{-5}' '10^{-4}' '10^{-3}'};
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('Na^+')

h(7) = subplot(5,2,7)  ; 
    [datlargest,datsmallest,newdat] = notails(crex.Fe2p_print,100);
    pcolor(Xr,ZZ,log10(crex.Fe2p_print));
    alcb = colorbar;
    alcb.Ticks = [ -10 -8 -6 -4];
    alcb.TickLabels = { '10^{-10}' '10^{-8}' '10^{-6}' '10^{-4}'};
    shading interp    
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('Fe^2^+')

h(8) = subplot(5,2,8)   ; 
    [datlargest,datsmallest,newdat] = notails(crex.Fe3p_print,100);
    pcolor(Xr,ZZ,log10(crex.Fe3p_print));
    alcb = colorbar;
    alcb.Ticks = [-12 -10 -8 -6 ];
    alcb.TickLabels = {'10^{-12}' '10^{-10}' '10^{-8}' '10^{-6}' };
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('Fe^3^+')
    xlabel('Day of Year')
    
 h(9) = subplot(5,2,9) ;
    [datlargest,datsmallest,newdat] = notails(crex.Al3p_print,100);
    pcolor(Xr,ZZ,log10(crex.Al3p_print));
    alcb = colorbar;
    alcb.Ticks = [-12 -10 -8 -6 -4];
    alcb.TickLabels = {'10^{-12}' '10^{-10}' '10^{-8}' '10^{-6}' '10^{-4}'};
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('Al^3^+')
    xlabel('Day of Year')

%% Percent change in concentration due to root exudation (Figure 7, 8, 9, (c) and Figure A.4)
figure()
h(1) = subplot(5,2,1)    ;
    pltdif = (crex.SiO2_print-crnoex.SiO2_print)./crnoex.SiO2_print.*100;
    pltdif(isnan(pltdif))=0;
    pltdif(isinf(pltdif))=0;
    pltdif = sign(pltdif).*log10(1+abs(pltdif)/10^-2);
    [datlargest,datsmallest,newdat] = notails(pltdif,100);
    customCMap1 = customcolormap(newdat,0,datlargest,datsmallest);
    pcolor(Xr,ZZ,pltdif);
    colormap(h(1),customCMap1)
    cbh = colorbar;
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('SiO_2_(_a_q_)')

h(2) = subplot(5,2,2) ;
    pltdif = (crex.TiO4H4_print-crnoex.TiO4H4_print)./crnoex.TiO4H4_print.*100;
    pltdif(isnan(pltdif))=0;
    pltdif(isinf(pltdif))=0;
     pltdif = sign(pltdif).*log10(1+abs(pltdif)/10^-4);
    [datlargest,datsmallest,newdat] = notails(pltdif,100);
    customCMap1 = customcolormap(newdat,0,datlargest,datsmallest);
    pcolor(Xr,ZZ,pltdif);
    colormap(h(2),customCMap1)
    cbh = colorbar;
    cbh.Ticks = [-2 -1 0 1];
    cbh.TickLabels = [{'-10^{-2}' '-10^{-3}' '0' '10^{-3}'}];
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('Ti(OH)_4_(_a_q_)')
    
h(3) = subplot(5,2,3)  ;
    pltdif = (crex.Ca2p_print-crnoex.Ca2p_print)./crnoex.Ca2p_print.*100;
    pltdif(isnan(pltdif))=0;
    pltdif(isinf(pltdif))=0;
     pltdif = sign(pltdif).*log10(1+abs(pltdif)/10^-2);
    [datlargest,datsmallest,newdat] = notails(pltdif,100);
    customCMap1 = customcolormap(newdat,0,datlargest,datsmallest);
    pcolor(Xr,ZZ,pltdif);
    colormap(h(3),customCMap1)
    cbh = colorbar;
    cbh.Ticks = [-1 0 1 2 3];
    cbh.TickLabels = num2cell([-0.1 0 0.1 1 10]);
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('Ca^2^+')

h(4) = subplot(5,2,4);
    pltdif = (crex.Mg2p_print-crnoex.Mg2p_print)./crnoex.Mg2p_print.*100;
    pltdif(isnan(pltdif))=0;
    pltdif(isinf(pltdif))=0;
     pltdif = sign(pltdif).*log10(1+abs(pltdif)/10^-2);
    [datlargest,datsmallest,newdat] = notails(pltdif,100);
    customCMap1 = customcolormap(newdat,0,datlargest,datsmallest);
    pcolor(Xr,ZZ,pltdif);
    colormap(h(4),customCMap1)
    cbh = colorbar;
    cbh.Ticks = [-2 0 2 4];
    cbh.TickLabels = num2cell([-1 0 1 100]);
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('Mg^2^+')

h(5) = subplot(5,2,5)  ;
    pltdif = (crex.Kp_print-crnoex.Kp_print)./crnoex.Kp_print.*100;
    pltdif(isnan(pltdif))=0;
    pltdif(isinf(pltdif))=0;
     pltdif = sign(pltdif).*log10(1+abs(pltdif)/10^-2);
    [datlargest,datsmallest,newdat] = notails(pltdif,100);
    customCMap1 = customcolormap(newdat,0,datlargest,datsmallest);
    pcolor(Xr,ZZ,pltdif);
    colormap(h(5),customCMap1)
    cbh = colorbar;
    cbh.Ticks = [-2 -1 0 1];
    cbh.TickLabels = num2cell([-1 -0.1 0 0.1]);
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('K^+')

h(6) = subplot(5,2,6)  ;
    pltdif = (crex.Nap_print-crnoex.Nap_print)./crnoex.Nap_print.*100;
    pltdif(isnan(pltdif))=0;
    pltdif(isinf(pltdif))=0;
     pltdif = sign(pltdif).*log10(1+abs(pltdif)/10^-8);
    [datlargest,datsmallest,newdat] = notails(pltdif,100);
    customCMap1 = customcolormap(newdat,0,datlargest,datsmallest);
    pcolor(Xr,ZZ,pltdif);
    colormap(h(6),customCMap1)
    cbh = colorbar;
    cbh.Ticks = [-2 0 2 4];
    cbh.TickLabels = [{'-10^{-6}' '0' '10^{-6}' '10^{-4}'}];
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('Na^+')

h(7) = subplot(5,2,7)  ; 
    pltdif = (crex.Fe2p_print-crnoex.Fe2p_print)./crnoex.Fe2p_print.*100;
    pltdif(isnan(pltdif))=0;
    pltdif(isinf(pltdif))=0;
     pltdif = sign(pltdif).*log10(1+abs(pltdif)/10^-2);
    [datlargest,datsmallest,newdat] = notails(pltdif,100);
    customCMap1 = customcolormap(newdat,0,datlargest,datsmallest);
    pcolor(Xr,ZZ,pltdif);
    colormap(h(7),customCMap1)
    cbh = colorbar;
    cbh.Ticks = [-2 0 2 4 ];
    cbh.TickLabels = num2cell([-1 0 1 100]);
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('Fe^2^+')

h(8) = subplot(5,2,8)   ; 
    pltdif = (crex.Fe3p_print-crnoex.Fe3p_print)./crnoex.Fe3p_print.*100;
    pltdif(isnan(pltdif))=0;
    pltdif(isinf(pltdif))=0;
     pltdif = sign(pltdif).*log10(1+abs(pltdif)/10^-4);
    [datlargest,datsmallest,newdat] = notails(pltdif,100);
    customCMap1 = customcolormap(newdat,0,datlargest,datsmallest);
    pcolor(Xr,ZZ,pltdif);
    colormap(h(8),customCMap1)
    cbh = colorbar;
    cbh.Ticks = [-5 -3 0 3 5 ];
    cbh.TickLabels = num2cell([-10 -0.1 0 0.1 10]);
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('Fe^3^+')
    xlabel('Day of Year')
    
 h(9) = subplot(5,2,9) ;
    pltdif = (crex.Al3p_print-crnoex.Al3p_print)./crnoex.Al3p_print.*100;
    pltdif(isnan(pltdif))=0;
    pltdif(isinf(pltdif))=0;
    pltdif = sign(pltdif).*log10(1+abs(pltdif)/10^-2);
    [datlargest,datsmallest,newdat] = notails(pltdif,100);
    customCMap1 = customcolormap(newdat,0,datlargest,datsmallest);
    pcolor(Xr,ZZ,pltdif);
    colormap(h(9),customCMap1)
    cbh = colorbar;
    cbh.Ticks = [-2 0 2 4];
    cbh.TickLabels = num2cell([-1 0 1 100]);
    shading interp
    set(gca, 'YDir','reverse')
    ax = gca;
    ax.Layer = 'top';
    ylabel('(m)')
    xlim([1 365])
    title('Al^3^+')
    xlabel('Day of Year')

set(findall(gcf,'-property','FontSize'),'FontSize',12)

%Log(concentration) vs Log(leaching) (Figure 7, 8, and 9 (b) and Figure A.5)
discharge = crex.qq_layer_day_point(13,:);
janmar = 1:90;
aprjun = 91:181;
julsep = 182:273;
octdec = 274:365;

figure()
subplot(6,2,1)
    loglog(discharge(janmar),crex.TLCH_SiO2_m2_print(janmar),'.','MarkerSize',12)
    hold on;
    loglog(discharge(aprjun),crex.TLCH_SiO2_m2_print(aprjun),'.','MarkerSize',12)
    hold on;
    loglog(discharge(julsep),crex.TLCH_SiO2_m2_print(julsep),'.','MarkerSize',12)
    hold on;
    loglog(discharge(octdec),crex.TLCH_SiO2_m2_print(octdec),'.','MarkerSize',12)
    title('SiO_2')
subplot(6,2,2)
    loglog(discharge(janmar),crex.TLCH_Nit_m2_print(janmar),'.','MarkerSize',12)
    hold on;
    loglog(discharge(aprjun),crex.TLCH_Nit_m2_print(aprjun),'.','MarkerSize',12)
    hold on;
    loglog(discharge(julsep),crex.TLCH_Nit_m2_print(julsep),'.','MarkerSize',12)
    hold on;
    loglog(discharge(octdec),crex.TLCH_Nit_m2_print(octdec),'.','MarkerSize',12)
    title('NO_3^-')
subplot(6,2,3)
    loglog(discharge(janmar),crex.TLCH_Al3p_m2_print(janmar),'.','MarkerSize',12)
    hold on;
    loglog(discharge(aprjun),crex.TLCH_Al3p_m2_print(aprjun),'.','MarkerSize',12)
    hold on;
    loglog(discharge(julsep),crex.TLCH_Al3p_m2_print(julsep),'.','MarkerSize',12)
    hold on;
    loglog(discharge(octdec),crex.TLCH_Al3p_m2_print(octdec),'.','MarkerSize',12)
    title('Al^3^+')
subplot(6,2,4)
    loglog(discharge(janmar),crex.TLCH_Amm_m2_print(janmar),'.','MarkerSize',12)
    hold on;
    loglog(discharge(aprjun),crex.TLCH_Amm_m2_print(aprjun),'.','MarkerSize',12)
    hold on;
    loglog(discharge(julsep),crex.TLCH_Amm_m2_print(julsep),'.','MarkerSize',12)
    hold on;
    loglog(discharge(octdec),crex.TLCH_Amm_m2_print(octdec),'.','MarkerSize',12)
    title('NH_4^+')
subplot(6,2,5)
    loglog(discharge(janmar),crex.TLCH_Ca2p_m2_print(janmar),'.','MarkerSize',12)
    hold on;
    loglog(discharge(aprjun),crex.TLCH_Ca2p_m2_print(aprjun),'.','MarkerSize',12)
    hold on;
    loglog(discharge(julsep),crex.TLCH_Ca2p_m2_print(julsep),'.','MarkerSize',12)
    hold on;
    loglog(discharge(octdec),crex.TLCH_Ca2p_m2_print(octdec),'.','MarkerSize',12)
    title('Ca^2^+')
subplot(6,2,6)
    loglog(discharge(janmar),crex.TLCH_Mg2p_m2_print(janmar),'.','MarkerSize',12)
    hold on;
    loglog(discharge(aprjun),crex.TLCH_Mg2p_m2_print(aprjun),'.','MarkerSize',12)
    hold on;
    loglog(discharge(julsep),crex.TLCH_Mg2p_m2_print(julsep),'.','MarkerSize',12)
    hold on;
    loglog(discharge(octdec),crex.TLCH_Mg2p_m2_print(octdec),'.','MarkerSize',12)
    title('Mg^2^+')
subplot(6,2,7)
    loglog(discharge(janmar),crex.TLCH_Kp_m2_print(janmar),'.','MarkerSize',12)
    hold on;
    loglog(discharge(aprjun),crex.TLCH_Kp_m2_print(aprjun),'.','MarkerSize',12)
    hold on;
    loglog(discharge(julsep),crex.TLCH_Kp_m2_print(julsep),'.','MarkerSize',12)
    hold on;
    loglog(discharge(octdec),crex.TLCH_Kp_m2_print(octdec),'.','MarkerSize',12)
    title('K^+')
subplot(6,2,8)
    loglog(discharge(janmar),crex.TLCH_Nap_m2_print(janmar),'.','MarkerSize',12)
    hold on;
    loglog(discharge(aprjun),crex.TLCH_Nap_m2_print(aprjun),'.','MarkerSize',12)
    hold on;
    loglog(discharge(julsep),crex.TLCH_Nap_m2_print(julsep),'.','MarkerSize',12)
    hold on;
    loglog(discharge(octdec),crex.TLCH_Nap_m2_print(octdec),'.','MarkerSize',12)
    title('Na^+')
subplot(6,2,9)
    loglog(discharge(janmar),crex.TLCH_Fe2p_m2_print(janmar),'.','MarkerSize',12)
    hold on;
    loglog(discharge(aprjun),crex.TLCH_Fe2p_m2_print(aprjun),'.','MarkerSize',12)
    hold on;
    loglog(discharge(julsep),crex.TLCH_Fe2p_m2_print(julsep),'.','MarkerSize',12)
    hold on;
    loglog(discharge(octdec),crex.TLCH_Fe2p_m2_print(octdec),'.','MarkerSize',12)
    title('Fe^2^+')
subplot(6,2,10)
    loglog(discharge(janmar),crex.TLCH_Fe3p_m2_print(janmar),'.','MarkerSize',12)
    hold on;
    loglog(discharge(aprjun),crex.TLCH_Fe3p_m2_print(aprjun),'.','MarkerSize',12)
    hold on;
    loglog(discharge(julsep),crex.TLCH_Fe3p_m2_print(julsep),'.','MarkerSize',12)
    hold on;
    loglog(discharge(octdec),crex.TLCH_Fe3p_m2_print(octdec),'.','MarkerSize',12)
    title('Fe^3^+')
h11 = subplot(6,2,11)
    loglog(discharge(janmar),crex.TLCH_TiO4H4_m2_print(janmar),'.','MarkerSize',12)
    hold on;
    loglog(discharge(aprjun),crex.TLCH_TiO4H4_m2_print(aprjun),'.','MarkerSize',12)
    hold on;
    loglog(discharge(julsep),crex.TLCH_TiO4H4_m2_print(julsep),'.','MarkerSize',12)
    hold on;
    loglog(discharge(octdec),crex.TLCH_TiO4H4_m2_print(octdec),'.','MarkerSize',12)
    title('Ti(OH)_4')
        legend([h11],{'Jan-Mar','Apr-Jun','Jul-Sep','Oct-Dec'},'Location','NorthEastOutside');







