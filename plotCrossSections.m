%Plots the emission and absorption cross sections

%clear all;
%close all;


%read in files and make sure wavelengths are in ascending order
%then resave file
%cs_absRAW = xlsread('silicate_abs.xlsx');
% cs_absRAW = xlsread('silica_abs_00.xlsx');
% [sorted_abs, sortIndex_abs] = unique(cs_absRAW(:,1));
% cs_absRAW = [sorted_abs cs_absRAW(sortIndex_abs,2)];
% %csvwrite('Las_Yb_06_01_abs.csv',cs_absRAW);
% csvwrite('silica_abs_00.csv',cs_absRAW);

%cs_emsRAW = xlsread('silicate_emi.xlsx');
% cs_emsRAW = xlsread('silica_emi_00.xlsx');
% [sorted_ems, sortIndex_ems] = unique(cs_emsRAW(:,1));
% cs_emsRAW = [sorted_ems cs_emsRAW(sortIndex_ems,2)];
% %csvwrite('Las_Yb_06_01_ems.csv',cs_emsRAW);
% csvwrite('silica_emi_00.csv',cs_emsRAW);

% wavelengths = (850:1150)*1e-9;
% cs_absRAW = xlsread('LAS_Yb_06_02_abs.xlsx');
% cs_abs = [wavelengths; interp1(cs_absRAW(:,1)*1e-9,cs_absRAW(:,2),wavelengths)].';
% cs_abs(isnan(cs_abs)) = 0;
% cs_emsRAW = xlsread('LAS_Yb_06_02_emi.xlsx');
% cs_ems = [wavelengths; interp1(cs_emsRAW(:,1)*1e-9,cs_emsRAW(:,2),wavelengths)].';
% cs_ems(isnan(cs_ems)) = 0;

%cross sectional areas 
wavelengths = (850:0.5:1100)*1e-9;
cs_absRAW = xlsread('abs_ZBLAN.xlsx');
cs_abs = [wavelengths; interp1(cs_absRAW(:,1),cs_absRAW(:,2),wavelengths)*1e-24].';
cs_abs(isnan(cs_abs)) = 0;
cs_emsRAW = xlsread('emm_ZBLAN.xlsx');
cs_ems = [wavelengths; interp1(cs_emsRAW(:,1),cs_emsRAW(:,2),wavelengths)*1e-24].';
cs_ems(isnan(cs_ems)) = 0;

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

%plot
figure(16)
hold on
plot(cs_abs(:,1)*1e9,cs_abs(:,2));
plot(cs_ems(:,1)*1e9,cs_ems(:,2));
xlabel('Wavelength (nm)');
ylabel('Cross Section (m^2)');
box on

%title('Cross section for Yb3+ doped Silica');
%legend('absorption','emission');

%calculate mean flourecent wavelength using Mina's formula
integral1 = cs_emsRAW(:,2)./cs_emsRAW(:,1).^4;%cs_emsRAW(:,1).*cs_emsRAW(:,2)./cs_emsRAW(:,1).^5;
integral2 = cs_emsRAW(:,2)./cs_emsRAW(:,1).^5;
lambdaF_SE = trapz(cs_emsRAW(:,1),integral1)/trapz(cs_emsRAW(:,1),integral2)*1e-9;
