clear all;
close all;

whatFun = 2;

global E11 E12 E13 E14 E21 E22 E23 E1 E2 g2 g1 tau21 N0 Area len core_radius loss;

%energy levels, m^-1
E11 = 0;
E12 = 206e2;
E13 = 290e2;
E14 = 380e2;
E21 = 10256e2;
E22 = 10417e2;
E23 = 10695e2;

E1 = [E11 E12 E13 E14];
E2 = [E21 E22 E23];

g2 = 3;
g1 = 4;

%Constants
tau21 = 1e-3; %s
N0 = 1.1e25 ; %m^-3;
core_radius = 1.5e-6; 
Area = pi*core_radius^2; %m^2, pump area
len = 10; %m; length of fiber
loss = 0.002;
h = 6.63e-34; %J*s
c = 3e8;
kT = 4.11e-21; %J

%cross sectional areas 
wavelengths = (870:1050)*1e-9;
cs_absRAW = xlsread('abs_ZBLAN.xlsx');
cs_abs = [wavelengths; interp1(cs_absRAW(:,1),cs_absRAW(:,2),wavelengths)*1e-24].';
cs_emsRAW = xlsread('emm_ZBLAN.xlsx');
cs_ems = [wavelengths; interp1(cs_emsRAW(:,1),cs_emsRAW(:,2),wavelengths)*1e-24].';

%calculate mean flourecent wavelength
integral1 = cs_emsRAW(:,1).*cs_emsRAW(:,2)./cs_emsRAW(:,1).^5;
integral2 = cs_emsRAW(:,2)./cs_emsRAW(:,1).^5;
lambda_f = trapz(cs_emsRAW(:,1),integral1)/trapz(cs_emsRAW(:,1),integral2);


%%% Analyze population inversion versus pump power %%%
if whatFun == 1
    %waveleghth of laser/cooling pump
    lambda_l = 980e-9;
    lambda_c = 1027e-9; 
    %cross sectional areas for a given wavelength
    indexl = find(round(cs_abs(:,1)*1e9 - lambda_l*1e9) == 0);
    indexc = find(round(cs_abs(:,1)*1e9 - lambda_c*1e9) == 0);
    cs_al = cs_abs(indexl,2);
    cs_el = cs_ems(indexl,2);
    cs_ac = cs_abs(indexc,2);
    cs_ec = cs_ems(indexc,2);
 
    %laser/cooling pump power
    P_l = 0.001:0.01:2; %W
    P_c = 0.00; %0:0.1:0.5; %W
    
    %does not take into account change in pump power across length of fiber
    %[N2,N1,N2j,N1i] = popInversionVsPower(lambda_l,lambda_c,cs_al,cs_el,cs_ac,cs_ec,P_l,P_c);
    %takes into accoutn change in pump power for ONE PUMP
    %[N2,N1,N2j,N1i] = popInversionVsPower2(lambda_l,cs_al,cs_el,P_l);
    %takes into account change in pump power for TWO PUMPS and calculates Q vs. P
    [N2,N1,N2j,N1i] = popInversionVsPower3(lambda_l,lambda_c,lambda_f,cs_al,cs_el,cs_ac,cs_ec,P_l,P_c);
end

%%% Analyze population inversion versus pump wavelength %%%
if whatFun == 2
    %laser/cooling pump power
    P_l = .0001; %W
    P_c = 0; %W
    
    %%does NOT account for change in power; TWO PUMPS
    %[N2,N1,N2j,N1i] = popInversionVsWavelength(P_l,P_c,cs_abs,cs_ems);
    %%accounts for change in power with ONE PUMP  and calculates Q vs. wavelength
    %[N2,N1,N2j,N1i,Q] = popInversionVsWavelength2(P_l,lambda_f,cs_abs,cs_ems);
    %%accounts for change in power with TWO PUMPS and calculates Q vs. wavelength
    [N2,N1,N2j,N1i] = popInversionVsWavelength3(P_l,P_c,lambda_f,cs_abs,cs_ems);
end

%calculate mean flourecence wavelength
sum1 = 4*N2j(1,150)+4*N2j(2,150)+4*N2j(3,150);
sum2 = N2j(1,150)*(E21-E11)+N2j(1,150)*(E21-E12)+N2j(1,150)*(E21-E13)+N2j(1,150)*(E21-E14)+...
    N2j(2,150)*(E22-E11)+N2j(2,150)*(E21-E12)+N2j(2,150)*(E21-E13)+N2j(2,150)*(E21-E14)+...
    N2j(3,150)*(E23-E11)+N2j(3,150)*(E23-E12)+N2j(3,150)*(E23-E13)+N2j(3,150)*(E23-E14);
lambda_f = sum1/sum2;






