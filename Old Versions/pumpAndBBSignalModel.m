%numerical result for change in power across fiber for one pump and broadboand signal

close all;
clear all;

len = 10;
dz = .5;
Pp0 = .2;
Psforward0 = 0;
Psbackward0 = [0.1209].'; %%GUESS THIS!!!!
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
dlam = 10e-9; %resolution of the signal spectrum
lambdaP = 935e-9;
lambdaS = [1013e-9:dlam:1013e-9]; 
freqP = c/lambdaP;
freqS = c./lambdaS;
dvS = c./(lambdaS-dlam/2)-c./(lambdaS+dlam/2);
%cross sectional areas for a given wavelength
indexP = find(round(cs_abs(:,1)*1e9 - lambdaP*1e9) == 0);
indexS = NaN(size(lambdaS));
for i = 1:length(lambdaS)
    indexS(i) = find(round(cs_abs(:,1)*1e9 - lambdaS(i)*1e9) == 0);
end
cs_aP = cs_abs(indexP,2);
cs_eP = cs_ems(indexP,2);
cs_aS = cs_abs(indexS,2).';
cs_eS = cs_ems(indexS,2).';

N0 = 1.1e25; %m^-3
f = 1/Area;
gamma = 1;
loss = 0.00;
IsatPa = h*freqP/cs_aP/tau21;
IsatSa = h*freqS./cs_aS/tau21;
IsatP = h*freqP/(cs_aP+cs_eP)/tau21;
IsatS = h*freqS./(cs_aS+cs_eS)/tau21;

z = 0:dz:len;  
Pp = zeros(1,length(z)); 
Psforward = zeros(length(lambdaS),length(z)); 
Psbackward = zeros(length(lambdaS),length(z)); 
Psbackward(:,1) = Psbackward0;
Pp(1) = Pp0; %W                                      % initial pump power


%N2(z) - Pstot is 1xlength(lambdaS) array where each element is the sum of
%the backward and forward power for that wavelength at position 
N2 = @(Pp,Pstot) (Pp*f/IsatPa+sum(Pstot*f./IsatSa))/(1+Pp*f/IsatP+sum(Pstot*f./IsatS))*N0;
%signal gain coefficients
gammaS = @(Pp,Pstot) N2(Pp,Pstot)*(cs_eS+cs_aS)-cs_aS*N0; %returns and array of size 1xlength(lambdaS)
gammaSE = @(Pp,Pstot) cs_eS*N2(Pp,Pstot); %returns and array of size 1xlength(lambdaS)
%pumps absorption coefficients
gammaP = @(Pp,Pstot) -1*N2(Pp,Pstot)*(cs_aP+cs_eP)+cs_aP*N0;

%pump power - Psf and Psb are 1xlength(lambdaS) arrays where each element
%is the forward and backward power (respectively) for that wavelegth and
%postion.
dPp_dz = @(z,Pp,Psf,Psb) -1*gammaP(Pp,Psf+Psb)*Pp;  
%forward signal power
dPsforward_dz = @(z,Pp,Psf,Psb) gammaS(Pp,Psf+Psb).*Psf+gammaSE(Pp,Psf+Psb)*2*h.*freqS.*dvS; %returns and array of size 1xlength(lambdaS)
%backward signal power
dPsbackward_dz = @(z,Pp,Psf,Psb) -1*gammaS(Pp,Psf+Psb).*Psb-gammaSE(Pp,Psf+Psb)*2*h.*freqS.*dvS; %returns and array of size 1xlength(lambdaS)

for i=1:(length(z)-1) % calculation loop
    
    RK1p = dPp_dz(z(i),Pp(i),Psforward(:,i).',Psbackward(:,i).');
    RK1f = dPsforward_dz(z(i),Pp(i),Psforward(:,i).',Psbackward(:,i).'); %array of 1xlength(lambdaS)
    RK1b = dPsbackward_dz(z(i),Pp(i),Psforward(:,i).',Psbackward(:,i).'); %array of 1xlength(lambdaS)
    
    RK2p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(:,i).'+0.5*dz*RK1f,Psbackward(:,i).'+0.5*dz*RK1b);
    RK2f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(:,i).'+0.5*dz*RK1f,Psbackward(:,i).'+0.5*dz*RK1b);
    RK2b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(:,i).'+0.5*dz*RK1f,Psbackward(:,i).'+0.5*dz*RK1b);
    
    RK3p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(:,i).'+0.5*dz*RK2f,Psbackward(:,i).'+0.5*dz*RK2b);
    RK3f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(:,i).'+0.5*dz*RK2f,Psbackward(:,i).'+0.5*dz*RK2b);
    RK3b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(:,i).'+0.5*dz*RK2f,Psbackward(:,i).'+0.5*dz*RK2b);
    
    RK4p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(:,i).'+0.5*dz*RK3f,Psbackward(:,i).'+0.5*dz*RK3b);
    RK4f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(:,i).'+0.5*dz*RK3f,Psbackward(:,i).'+0.5*dz*RK3b);
    RK4b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(:,i).'+0.5*dz*RK3f,Psbackward(:,i).'+0.5*dz*RK3b);
    
    Pp(i+1) = Pp(i) + (1/6)*(RK1p+2*RK2p+2*RK3p+RK4p)*dz;
    Psforward(:,i+1) = Psforward(:,i) + (1/6)*(RK1f+2*RK2f+2*RK3f+RK4f).'*dz;  % main equation
    Psbackward(:,i+1) = Psbackward(:,i) + (1/6)*(RK1b+2*RK2b+2*RK3b+RK4b).'*dz;  % main equation
end



%graphs
figure(1)
grid on
plot(z,Pp);
xlabel('z (m)');
ylabel('Pp(z) (W)');
title('Change in pump Power along fiber');

figure(2)
hold on
grid on
for i = 1:length(lambdaS)
    plot(z,Psforward(i,:),'DisplayName',sprintf('%g nm',lambdaS(i)*1e9));
end
xlabel('z (m)');
ylabel('Psforward(z) (W)');
title('Change in forward signal along fiber');
legend('show');

figure(3)
hold on
grid on
for i = 1:length(lambdaS)
    plot(z,Psbackward(i,:),'DisplayName',sprintf('%g nm',lambdaS(i)*1e9));
end
xlabel('z (m)');
ylabel('Psbackward(z) (W)');
title('Change in backward signal along fiber');
legend('show')

figure(4)
hold on 
grid on
plot(z,Pp);
plot(z,sum(Psforward));
plot(z,sum(Psbackward));
xlabel('z (m)');
ylabel('Power (W)');
title('Change in Power along fiber');
legend('pump','forward signal','backward signal');

figure(5)
hold on
grid on
PsforwardNorm = Psforward(:,length(lambdaS))/max(Psforward(:,length(lambdaS)));
PsbackwardNorm = Psbackward(:,length(lambdaS))/max(Psbackward(:,length(lambdaS)));
plot(lambdaS*1e9,PsforwardNorm);
plot(lambdaS*1e9,PsbackwardNorm);
xlabel('wavelength (nm)');
ylabel('power (nm)');
title('Normalized Spectrum of forward and backward signal')
legend('forward signal','backward signal');






