% takes in an initial temperature distribution and iteratively solves for the SS temp distribution
% taking into account convection and conduction.

%9/18/17 - probably doesn't work. Tried to assess functionality with
%'tempDistributionTestFin.m' to see if it would converge to know solution
%and it doesn't...so this code probably doesn't work either....

%close all;

Lyb = len;
dzYb = dz;

maxIterations = 500;
convError = 0.00001; %end when interations change by less than 1%

h = 81.4; %W/m^2/k
k = 1.4e3; %W/k
r_core = 3.1;%e-6; %m
r_clad = 67.3e-6; %m
peri = 2*pi*r_clad; %perimeter of fiber
Ac = pi*r_clad^2; %cross sectional area of fiber
m = sqrt(h*peri/k/Ac);
RT = 300; %K - Room Temperature
L1 = 20; %m - length of SMF fiber at beginning
dz1 = L1/100;
L2 = 20; %m - length of SMF fiber at end
dz2 = L2/100;

z = [0:dz1:L1 L1:dzYb:(L1+Lyb) (L1+Lyb+dz2):dz2:(L1+Lyb+L2)];
Ybtemp = RT+dTz;
initialT1 = RT+dTz(1)*exp(-m*(L1-z(1:floor(L1/dz1))));
initialTYb = Ybtemp;
initialT2 = RT+dTz(end)*exp(-m*(z(end-floor(L2/dz2):end)-L1-Lyb));
initialTemp = [initialT1 initialTYb initialT2];

nextTemp = zeros(size(initialTemp));

iteration = 1;

while iteration <= maxIterations
 nextTemp(1) = initialTemp(1);
 ai = 2;
 bi = 1;
 ci = 1;
 di = [dz1^2/k*h*(initialT1-RT) dzYb^2/k*(h*(initialTYb-RT)+dQz_dt2) dz2^2/k*h*(initialT2-RT)];
 %di = zeros(1,302);
 for i = 2:length(z)-1
     nextTemp(i) = (bi*initialTemp(i+1)+ci*initialTemp(i-1)+di(i))/ai;
 end
 nextTemp(end) = initialTemp(end);
 
 %%check boundary conditions
 if sum(abs(initialTemp-nextTemp)./initialTemp > convError) == 0 
     break
 end
 
 figure(1)
 hold on
 plot(z,initialTemp);
 xlabel('z (m)');
 ylabel('temperature (K)');
 
 initialTemp = nextTemp;
 initialT1 = nextTemp(1:length(initialT1));
 initialTYb = nextTemp(length(initialT1)+1:length(initialT1)+length(initialTYb));
 initialT2 = nextTemp(end-length(initialT2)+1:end);
 
 iteration = iteration + 1;

end

figure(2)
hold on
plot(z,nextTemp);
xlabel('z (m)');
ylabel('temperature (K)');
title('SS temperature distribution');


dQ = [ dQz_dt2 zeros(1,L2/dz2)];
