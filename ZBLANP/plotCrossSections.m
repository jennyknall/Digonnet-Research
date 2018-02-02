%Plots the emission and absorption cross sections

clear all;
close all;


%read in files and make sure wavelengths are in ascending order
%then resave file
% cs_absRAW = xlsread('emm_ZBLANPLei.xlsx');
% [sorted_abs, sortIndex_abs] = unique(cs_absRAW(:,1));
% cs_absRAW = [sorted_abs cs_absRAW(sortIndex_abs,2)];
% cs_absRAW(:,1) = cs_absRAW(:,1)*1e-9;
% xlswrite('emm_ZBLANPLei',cs_absRAW);


%cs_emsRAW = xlsread('silicate_emi.xlsx');
% cs_emsRAW = xlsread('silica_emi_00.xlsx');
% [sorted_ems, sortIndex_ems] = unique(cs_emsRAW(:,1));
% cs_emsRAW = [sorted_ems cs_emsRAW(sortIndex_ems,2)];
% %csvwrite('Las_Yb_06_01_ems.csv',cs_emsRAW);
% csvwrite('silica_emi_00.csv',cs_emsRAW);

%cross sectional areas - Silica
wavelengths = (850:1150)*1e-9;
cs_absRAWSilica = xlsread('LAS_Yb_06_02_abs.xlsx');
cs_absSilica = [wavelengths; interp1(cs_absRAWSilica(:,1)*1e-9,cs_absRAWSilica(:,2),wavelengths)].';
cs_absSilica(isnan(cs_absSilica)) = 0;
cs_emsRAWSilica = xlsread('LAS_Yb_06_02_emi.xlsx');
cs_emsSilica = [wavelengths; interp1(cs_emsRAWSilica(:,1)*1e-9,cs_emsRAWSilica(:,2),wavelengths)].';
cs_emsSilica(isnan(cs_emsSilica)) = 0;
cs_emsRAW_MCSilica = xlsread('LAS_Yb_06_02_emi.xlsx');
cs_ems_MCSilica = [wavelengths; interp1(cs_emsRAW_MCSilica(:,1),cs_emsRAW_MCSilica(:,2),wavelengths)].';
cs_ems_MCZSilica(isnan(cs_ems_MCSilica)) = 0;

%cross sectional areas - ZBLANP
wavelengths = (850:1:1100)*1e-9;
cs_absRAWZBLANP = xlsread('abs_ZBLANP_MC.xlsx');
cs_absZBLANP = [wavelengths; interp1(cs_absRAWZBLANP(:,1),cs_absRAWZBLANP(:,2),wavelengths)].';
cs_absZBLANP(isnan(cs_absZBLANP)) = 0;
cs_emsRAWZBLANP = xlsread('emm_ZBLANP_MC.xlsx');
cs_emsZBLANP = [wavelengths; interp1(cs_emsRAWZBLANP(:,1),cs_emsRAWZBLANP(:,2),wavelengths)].';
cs_emsZBLANP(isnan(cs_emsZBLANP)) = 0;
cs_emsRAW_MCZBLANP = xlsread('emm_ZBLANP_MC.xlsx');
cs_ems_MCZBLANP = [wavelengths; interp1(cs_emsRAW_MCZBLANP(:,1),cs_emsRAW_MCZBLANP(:,2),wavelengths)].';
cs_ems_MCZBLANP(isnan(cs_ems_MCZBLANP)) = 0;

% csvwrite('abs_ZBLANPLei_Digi.csv',cs_absZBLANP);
% csvwrite('emm_ZBLANPLei_Digi.csv',cs_emsZBLANP);


%cross sectional areas - ZBLAN
wavelengths = (850:1.5:1100)*1e-9;
cs_absRAWZBLAN = xlsread('abs_ZBLAN.xlsx');
cs_absZBLAN = [wavelengths; interp1(cs_absRAWZBLAN(:,1),cs_absRAWZBLAN(:,2),wavelengths)*1e-24].';
cs_absZBLAN(isnan(cs_absZBLAN)) = 0;
cs_emsRAWZBLAN = xlsread('emm_ZBLAN.xlsx');
cs_emsZBLAN = [wavelengths; interp1(cs_emsRAWZBLAN(:,1),cs_emsRAWZBLAN(:,2),wavelengths)*1e-24].';
cs_emsZBLAN(isnan(cs_emsZBLAN)) = 0;
cs_emsRAW_MCZBLAN = xlsread('emmMC_ZBLAN.xlsx');
cs_ems_MCZBLAN = [wavelengths; interp1(cs_emsRAW_MCZBLAN(:,1),cs_emsRAW_MCZBLAN(:,2),wavelengths)].';
cs_ems_MCZBLAN(isnan(cs_ems_MCZBLAN)) = 0;

% numWaves =length(wavelengths);
% lambda_avg = 0;
% for i = 1:numWaves
%     lambda_avg = lambda_avg+wavelengths(i)*cs_ems(i,2);
% end
% lambda_avg = lambda_avg/sum(cs_ems(:,2))



%cross sectional areas 
%wavelengths = (850:1090)*1e-9;
% wavelengths = (850:1080)*1e-9;
% cs_abs = [wavelengths; interp1(cs_absRAW(:,1)*1e-9,cs_absRAW(:,2),wavelengths)].';
% cs_abs(isnan(cs_abs)) = 0;
% % cs_abs = [cs_abs(1:end-5,1) cs_abs(1:end-5,2)];   
% %cs_abs(:,2) = smooth(cs_abs(:,2)); 
% cs_ems = [wavelengths; interp1(cs_emsRAW(:,1)*1e-9,cs_emsRAW(:,2),wavelengths)].';
% cs_ems(isnan(cs_ems)) = 0;
% emission = cs_ems(:,2);
% newEmission = [emission(1:150); smooth(emission(151:235),15)];
% cs_ems = [cs_ems(1:end-5,1) newEmission(1:end-5)];   

% 
% csvwrite('silicate_abs.csv',cs_abs);
% csvwrite('silicate_emi.csv',cs_ems);

%plot - ZBLAN
figure(1)
hold on
plot(cs_absZBLAN(:,1)*1e9,cs_absZBLAN(:,2));
plot(cs_emsZBLAN(:,1)*1e9,cs_emsZBLAN(:,2));
plot(cs_ems_MCZBLAN(:,1)*1e9,cs_ems_MCZBLAN(:,2));
xlabel('Wavelength (nm)');
ylabel('Cross Section (m^2)');
title('ZBLAN')
legend('absorption cross section','given emission cross section','McCumber emission');
box on

%plot - ZBLANP
figure(2)
hold on
plot(cs_absZBLANP(:,1)*1e9,cs_absZBLANP(:,2)*1e24);
plot(cs_emsZBLANP(:,1)*1e9,cs_emsZBLANP(:,2)*1e24);
%plot(cs_ems_MCZBLANP(:,1)*1e9,cs_ems_MCZBLANP(:,2));
xlabel('Wavelength (nm)');
ylabel('Cross Section (10^{-24} m^2)');
title('ZBLANP')
%legend('absorption cross section','given emission cross section');%,'McCumber emission');
box on
xlim([870 1100])
ylim([0 1.4]);

%plot - Silica
figure(3)
hold on
plot(cs_absSilica(:,1)*1e9,cs_absSilica(:,2));
plot(cs_emsSilica(:,1)*1e9,cs_emsSilica(:,2));
plot(cs_ems_MCZSilica(:,1)*1e9,cs_ems_MCSilica(:,2));
xlabel('Wavelength (nm)');
ylabel('Cross Section (m^2)');
title('Silica')
legend('absorption cross section','given emission cross section','McCumber emission');
box on

%title('Cross section for Yb3+ doped Silica');
%legend('absorption','emission');

% %calculate mean flourecent wavelength using Mina's formula
% integral1 = cs_emsRAW(:,2)./cs_emsRAW(:,1).^4;%cs_emsRAW(:,1).*cs_emsRAW(:,2)./cs_emsRAW(:,1).^5;
% integral2 = cs_emsRAW(:,2)./cs_emsRAW(:,1).^5;
% lambdaF_SE = trapz(cs_emsRAW(:,1),integral1)/trapz(cs_emsRAW(:,1),integral2)*1e-9;
