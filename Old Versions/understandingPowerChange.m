% numerical and analytical solutions to dp/dz for one pump

%close all;
clear all;

len = 5;
dz = 0.5;
y0 = 0.005;
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
lambda = 1010e-9;
freq = c/lambda;
%cross sectional areas for a given wavelength
index = find(round(cs_abs(:,1)*1e9 - lambda*1e9) == 0);
cs_a = cs_abs(index,2);
cs_e = cs_ems(index,2);

N0 = 1.1e25 ; %m^-3
f = 1/Area;
gamma = 1;
loss = 0.00;
Isat = h*freq/(cs_a+cs_e)/tau21;

A = (cs_a+cs_e)*f/h/freq*cs_a*tau21*gamma*N0;
B = f/Isat;
C = cs_a*N0*gamma+loss*gamma;

%%numerical solution
z = 0:dz:len;                                         
Ynum = zeros(1,length(z)); 
Ynum(1) = y0; %W                                      % initial pump power
dy_dz = @(z,y) ((A-C*B)*y^2-C*y)/(1+B*y);  

for i=1:(length(z)-1)                              % calculation loop
    k_1 = dy_dz(z(i),Ynum(i));
    k_2 = dy_dz(z(i)+0.5*dz,Ynum(i)+0.5*dz*k_1);
    k_3 = dy_dz((z(i)+0.5*dz),(Ynum(i)+0.5*dz*k_2));
    k_4 = dy_dz((z(i)+dz),(Ynum(i)+k_3*dz));

    Ynum(i+1) = Ynum(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dz;  % main equation
end

%%analytical solution
%if loss DOES equal 0, run this code (since A-C*D == 0):
const = 1/C*log(y0)+B/C*y0;
fun = @(y,z) -1/C*log(y)-B/C*y+const-z;
Yana = [];
for z = 0:dz:len
    myfun = @(y) fun(y,z);
    Yana = [Yana fzero(myfun,[1e-200 y0])];
end

% %if loss does NOT equal 0, run this code (since A-C*D ~= 0):
% const = -B/(A-C*B)*log(y0)-(1+B*C/(A-C*B))/C*(log(abs(C-(A-C*B)*y0))-log(y0));
% fun = @(y,z) B/(A-C*B)*log(y)+(1+B*C/(A-C*B))/C*(log(abs(C-(A-C*B)*y))-log(y))+const-z;
% Yana = [];
% for z = 0:dz:len
%     myfun = @(y) fun(y,z);
%     Yana = [Yana fzero(myfun,[1e-200 y0])];
% end

%graphs
figure(1)
z = 0:dz:len;
hold on
plot(z,Ynum);
plot(z,Yana);
%semilogy(z,Ynum,z,Yana);
%plot(z,Ynum)
grid on
xlabel('z (m)');
ylabel('P(z) (W)');
title(sprintf('Change in Power along fiber for P0 = %g W',y0));
legend('numerical solution','analytical solution');


