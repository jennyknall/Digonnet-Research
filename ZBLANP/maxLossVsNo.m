close all
clear all

N0_wtPercent = [0.00005:.00005:.001 .001:.005:.1 .1:.01:4];
N0 = N0_wtPercent*2.42e26;  %m^-3

core_radius = 3.1e-6; 
tauRad = 1.7e-3;
Nc = 6.47e27; %m^-3

%cross sectional areas 
wavelengths = (870:1100)*1e-9;
cs_absRAW = xlsread('abs_ZBLANP_MC.xlsx');
cs_abs = [wavelengths; interp1(cs_absRAW(:,1),cs_absRAW(:,2),wavelengths)].';
cs_abs(isnan(cs_abs)) = 0;
cs_emsRAW = xlsread('emm_ZBLANP_MC.xlsx');
cs_ems = [wavelengths; interp1(cs_emsRAW(:,1),cs_emsRAW(:,2),wavelengths)].';
cs_ems(isnan(cs_ems)) = 0;

%calculate mean flourecent wavelength using Mina's formula
integral1 = cs_emsRAW(:,2)./cs_emsRAW(:,1).^4;%cs_emsRAW(:,1).*cs_emsRAW(:,2)./cs_emsRAW(:,1).^5;
integral2 = cs_emsRAW(:,2)./cs_emsRAW(:,1).^5;
lambdaF_SE = trapz(cs_emsRAW(:,1),integral1)/trapz(cs_emsRAW(:,1),integral2);

%waveleghth of laser/cooling pump
lambdaP = 1015e-9;
%cross sectional areas for a given wavelength
indexP = find(round(cs_abs(:,1)*1e9 - lambdaP*1e9) == 0);
cs_aP = cs_abs(indexP,2);
cs_eP = cs_ems(indexP,2);

lam = 1000e-9;
NA = .13; 
V = 2*pi*core_radius/lam*NA;
w = core_radius*(0.65+1.619/V^1.5+2.879/V^6);
etta = 1-exp(-2*core_radius^2/w^2);


tauNonRad = @(N) 2*pi/9*tauRad*(Nc./N).^2;
tau = @(N) tauRad*tauNonRad(N)./(tauRad+tauNonRad(N));

maxLoss = 10/log(10)*etta*cs_aP*N0.*(tau(N0)/tauRad*lambdaP/lambdaF_SE-1);

figure(1)
hold on
box on
plot(N0_wtPercent,maxLoss)
xlabel('Yb concentration (wt% Yb)');
ylabel('Max absorptive loss (dB/m)');
ylim([0 max(maxLoss)])





