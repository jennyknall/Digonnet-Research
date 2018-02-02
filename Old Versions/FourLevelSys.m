clear all; close all;

%intrinsic constants
h = 6.63e-34;
c = 3e8;
E0 = 0;
E1 = .5e5; %m^-1
E2 = 9.5e5; %m^-1
E3 = 10e5; %m^-1
lam23 = 1/(E3-E2); %m
lam01 = 1/(E1-E0); %m
lam12 = 1/(E2-E1); %m
lam03 = 1/(E3-E0); %m
delEu = h*c/lam23;
delEg = h*c/lam01;
wu = 1e12;
wg = 1e12;
tr = 1e-3;
kT = 4.11e-21; %J
pumpCS = 10e-24;

%Controlable constants
Area = 1e-10; %m^2
vp = c/lam12;
P = 0.001:0.001:1;
Wp = pumpCS/h/vp*P/Area;
Ntot = 2e26; %m^-3

A = (Wp+(1/tr-wg*exp(-delEg/kT))*(1+wu*exp(-delEu/kT)/(wu+2/tr)))./(Wp+wg+wg*exp(-delEg/kT));
B = wu*exp(-delEu/kT)/(wu+2/tr);
C = Ntot*wg*exp(-delEg/kT)./(Wp+wg+wg*exp(-delEg/kT));

N2 = Wp.*C./(-Wp.*A+wu*exp(-delEu/kT)-B*wu+2/tr+Wp);
N3 = N2*B;
N1 = N2.*A+C;
N0 = Ntot-N1-N2-N3;

figure(1)
hold on
%plot(P,N0)
plot(P,N1)
plot(P,N2)
plot(P,N3)
plot(P,N3./N2)
%legend('N0','N1','N2','N3','N3/N2')
legend('N1','N2','N3','N3/N2')
title('Population density vs. Power')
xlabel('Power (W)')
ylabel('Population Density (m^-3)')

figure(2)
plot(P,N2./N1)
title('Population inversion vs. Power')
xlabel('Power (W)')
ylabel('N2/N1')

%%Flourescent Intensity/time

Iout = N3/tr*(h*c*(E3-E0)+h*c*(E3-E1))+N2/tr*(h*c*(E2-E0)+h*c*(E2-E1));
Iin = pumpCS*P/Area.*(N1-N2);

delI = Iin-Iout;
