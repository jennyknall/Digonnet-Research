%numerical result for change in power across fiber for two pumps

%close all;
clear all;

len = 5;
dz = .5;
y0 = .005;
x0 = 0.000;
core_radius = 1.5e-6; 
Area = pi*core_radius^2; %m^2, pump area
tau21 = 1e-3;
c = 3e8;
h = 6.63e-34; %J*s

% %cross sectional areas 
wavelengths = (870:1050)*1e-9;
cs_absRAW = xlsread('abs_ZBLAN.xlsx');
cs_abs = [wavelengths; interp1(cs_absRAW(:,1),cs_absRAW(:,2),wavelengths)*1e-24].';
cs_emsRAW = xlsread('emm_ZBLAN.xlsx');
cs_ems = [wavelengths; interp1(cs_emsRAW(:,1),cs_emsRAW(:,2),wavelengths)*1e-24].';
% 
% %waveleghth of laser/cooling pump
lambday = 1010e-9;
lambdax = 1027e-9; 
freqy = c/lambday;
freqx = c/lambdax;
%cross sectional areas for a given wavelength
indexy = find(round(cs_abs(:,1)*1e9 - lambday*1e9) == 0);
indexx = find(round(cs_abs(:,1)*1e9 - lambdax*1e9) == 0);
cs_ay = cs_abs(indexy,2);
cs_ey = cs_ems(indexy,2);
cs_ax = cs_abs(indexx,2);
cs_ex = cs_ems(indexx,2);

N0 = 1.1e25; %m^-3
f = 1/Area;
gamma = 1;
loss = 0.00;
Isaty = h*freqy/(cs_ay+cs_ey)/tau21;
Isatx = h*freqx/(cs_ax+cs_ex)/tau21;

Aprime = cs_ay+cs_ey;
Bprime = f*cs_ay/h/freqy*tau21*N0*gamma;
Cprime = f*cs_ax/h/freqx*tau21*N0*gamma;
Dprime = f/Isaty;
Eprime = f/Isatx;
Fprime = cs_ay*N0*gamma+loss*gamma;
Gprime = cs_ax+cs_ex;
Hprime = cs_ax*N0*gamma+loss*gamma;

A = Aprime*Bprime-Fprime*Dprime;
B = Aprime*Cprime-Fprime*Eprime;
C = Fprime;
F = Gprime*Cprime-Hprime*Eprime;
G = Gprime*Bprime-Hprime*Dprime;
H = Hprime;
D = Dprime;
E = Eprime;

z = 0:dz:len;                                         
Y = zeros(1,length(z)); 
X = zeros(1,length(z));
Y(1) = y0; %W                                      % initial laser pump power
X(1) = x0; %W                                      % initial cooling pump power
%laser pumping power
dY_dz = @(z,y,x) (A*y^2+B*x*y-C*y)/(D*y+E*x+1);  
%cooling pumping power
dX_dz = @(z,y,x) (F*x^2+G*x*y-H*x)/(D*y+E*x+1);

for i=1:(length(z)-1)                              % calculation loop
    k_1 = dY_dz(z(i),Y(i),X(i));
    j_1 = dX_dz(z(i),Y(i),X(i));
    k_2 = dY_dz(z(i)+0.5*dz,Y(i)+0.5*dz*k_1,X(i)+0.5*dz*j_1);
    j_2 = dX_dz(z(i)+0.5*dz,Y(i)+0.5*dz*k_1,X(i)+0.5*dz*j_1);
    k_3 = dY_dz((z(i)+0.5*dz),(Y(i)+0.5*dz*k_2),(X(i)+0.5*dz*j_2));
    j_3 = dX_dz((z(i)+0.5*dz),(Y(i)+0.5*dz*k_2),(X(i)+0.5*dz*j_2));
    k_4 = dY_dz((z(i)+dz),(Y(i)+k_3*dz),(X(i)+j_3*dz));
    j_4 = dX_dz((z(i)+dz),(Y(i)+k_3*dz),(X(i)+j_3*dz));

    Y(i+1) = Y(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dz;  % main equation
    X(i+1) = X(i) + (1/6)*(j_1+2*j_2+2*j_3+j_4)*dz;  % main equation
end

%alternate way of calculating ODE solutions using matlabs built in RK
ODEs = @(z,p) [dY_dz(z,p(1),p(2));dX_dz(z,p(1),p(2))];
[zee, sols] = ode45(ODEs,[0 len],[y0 x0]);

%graphs
figure(3)
plot(z,Y);
xlabel('z (m)');
ylabel('Pl(z) (W)');
title('Change in laser pump Power along fiber');

figure(4)
plot(z,X);
xlabel('z (m)');
ylabel('Pc(z) (W)');
title('Change in cooling pump Power along fiber');


