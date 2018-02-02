% takes in an initial temperature distribution and iteratively solves for the SS temp distribution
% taking into account convection and conduction.

% of a cylindical fin

%9/18/17 - code does NOT converge even when the analytical solution is put
%as the inital guess. For smaller dz, the code thinks that it has converged
%but it doesn't really change much from the original guess...

close all;

L = 1;
dz = L/500;
z = 0:dz:L;

maxIterations = 500;
convError = 0.00001; %end when interations change by less than 1%

h = 81.4; %W/m^2/k
k = 1.4e3; %W/k
r = 1e-2;
peri = 2*pi*r; %perimeter of fiber
Ac = pi*r^2; %cross sectional area of fiber
m = sqrt(h*peri/k/Ac);
T1 = 350;
RT = 300; %K - Room Temperature

initialTemp = -z*(T1-RT)/L+T1;

nextTemp = zeros(size(initialTemp));

iteration = 1;

while iteration <= maxIterations
 nextTemp(1) = T1;
 ai = 2;
 bi = 1;
 ci = 1;
 di = dz^2/k*h*(initialTemp-RT);
 %di = zeros(1,302);
 for i = 2:length(z)-1
     %di(i) = dz^2/k*h*(nextTemp(i-1)-RT);
     nextTemp(i) = (bi*initialTemp(i+1)+ci*nextTemp(i-1)+di(i))/ai;
 end
 aN = 1;
 cN = 1/(1+h*dz/k);
 dN = dz/k*h*RT/(1+h*dz/k);
 nextTemp(end) = (cN*initialTemp(end-1)+dN)/aN;
 
 
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
 
 iteration = iteration + 1;

end

solTemp = (cosh(m*L-z)+h/m/k*sinh(m*(L-z)))/(cosh(m*L)+h/m/k*sinh(m*L))*(T1-RT)+RT;

figure(2)
hold on
plot(z,nextTemp);
plot(z,solTemp);
xlabel('z (m)');
ylabel('temperature (K)');
title('SS temperature distribution');
legend('numerical sol','analytical sol');

