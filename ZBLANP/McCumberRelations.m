close all

%cross sectional areas - ZBLAN
wavelengths = (850:1:1100)*1e-9;
cs_absRAWZBLAN = xlsread('abs_ZBLANPLei_Digi.xlsx');
%cs_absRAWZBLAN = xlsread('LAS_Yb_06_02_abs.xlsx');
%cs_absZBLAN = [wavelengths; interp1(cs_absRAWZBLAN(:,1)*1e-9,cs_absRAWZBLAN(:,2),wavelengths)].'; %for silica cs
cs_absZBLAN = [wavelengths; interp1(cs_absRAWZBLAN(:,1),cs_absRAWZBLAN(:,2),wavelengths)].'; %for ZBLANP cs
cs_absZBLAN(isnan(cs_absZBLAN)) = 0;
cs_emsRAWZBLAN = xlsread('emm_ZBLANPLei_Digi.xlsx');
%cs_emsRAWZBLAN = xlsread('LAS_Yb_06_02_emi.xlsx');
%cs_emsZBLAN = [wavelengths; interp1(cs_emsRAWZBLAN(:,1)*1e-9,cs_emsRAWZBLAN(:,2),wavelengths)].'; %for silica cs
cs_emsZBLAN = [wavelengths; interp1(cs_emsRAWZBLAN(:,1),cs_emsRAWZBLAN(:,2),wavelengths)].'; %for ZBLANP cs
cs_emsZBLAN(isnan(cs_emsZBLAN)) = 0;

e2 = [10254 10425 10668]*1e2; %1/m
e2m = [10254*1.98630e-23 10425*1.98630e-23 10668*1.98630e-23];
e1l = [0,181*1.98630e-23,321*1.98630e-23,482*1.98630e-23]; 
%Subscript[e, 2m]=List[10256*1.98630*10^-23,10417*1.98630*10^-23,10695*1.98630*10^-23];
%Subscript[e, 1l]=List[0 206*1.98630e-23 290*1.98630e-23 380*1.98630e-23]; *) 

c = 3e8;
K=1.3807e-23;
T= 273;
h = 6.63e-34;

e0 = -(e2m(1)-e1l(1));
Z1 = 0;
Z2 = 0;
for i = 1:4
    Z1 = Z1 + exp(-(e1l(i)-e1l(1))/(K*T));
end 
for i = 1:3
    Z2 = Z2 + exp(-(e2m(i)-e2m(1))/(K*T));
end

%epsilon = -log(Z1/Z2*exp(e0/(K*T)))*(K*T);
epsilon = h*c/974e-9;%e2(1)

emm_MC= cs_absZBLAN(:,2).*exp((epsilon-h*c./cs_absZBLAN(:,1))/(K*T));
abs_MC= cs_emsZBLAN(:,2).*exp((h*c./cs_emsZBLAN(:,1)-epsilon)/(K*T));

figure(1)
hold on
plot(cs_absZBLAN(:,1)*1e9,cs_absZBLAN(:,2)*1e24);
plot(cs_absZBLAN(:,1)*1e9,abs_MC*1e24);
xlabel('Wavelength (nm)')
ylabel('Cross sections (10^{-24} m^2)')
title('absorption')
legend('measured','McCumber')
xlim([870 1100]);
ylim([0 1.4]);
box on 

figure(2)
hold on
plot(cs_emsZBLAN(:,1)*1e9,cs_emsZBLAN(:,2)*1e24);
plot(cs_emsZBLAN(:,1)*1e9,emm_MC*1e24);
xlabel('Wavelength (nm)')
ylabel('Cross sections (10^{-24} m^2)')
title('emission')
legend('measured','McCumber')
xlim([870 1100]);
ylim([0 1.4]);
box on 

newEmm_MC = [cs_emsZBLAN(:,1) [emm_MC(1:126); cs_emsZBLAN(127:end,2)]];
newAbs_MC = [cs_absZBLAN(:,1) [cs_absZBLAN(1:125,2); abs_MC(126:end)]];

figure(3)
hold on
plot(newEmm_MC(:,1)*1e9,newEmm_MC(:,2)*1e24);
plot(cs_emsZBLAN(:,1)*1e9,cs_emsZBLAN(:,2)*1e24);
plot(newAbs_MC(:,1)*1e9,newAbs_MC(:,2)*1e24);
plot(cs_absZBLAN(:,1)*1e9,cs_absZBLAN(:,2)*1e24);
xlabel('Wavelength (nm)')
ylabel('Cross sections (10^{-24} m^2)')
title('McCumber Relations')
legend('MC emm','Meas. emm','MC abs','Meas. abs')
xlim([870 1100]);
ylim([0 1.4]);
box on 

% xlswrite('emm_ZBLANP_MC_peak',newEmm_MC);
% xlswrite('abs_ZBLANP_MC_peak',newAbs_MC);





