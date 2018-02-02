%numerical result for change in power across fiber for one pump and one
%signal. Automatically converges. 

close all;
clear all;

len = 5;
dz = .1;
Pp0 = 0.2;
Psforward0 = 0;
Psbackward0 = .2; %%GUESS THIS!!!!
maxError = 1e-7;
maxIterations = 30;

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
lambdaP = 935e-9;
lambdaS = 1010e-9; 
freqP = c/lambdaP;
freqS = c/lambdaS;
dvS = c/(lambdaS-.5e-9)-c/(lambdaS+.5e-9);
%cross sectional areas for a given wavelength
indexP = find(round(cs_abs(:,1)*1e9 - lambdaP*1e9) == 0);
indexS = find(round(cs_abs(:,1)*1e9 - lambdaS*1e9) == 0);
cs_aP = cs_abs(indexP,2);
cs_eP = cs_ems(indexP,2);
cs_aS = cs_abs(indexS,2);
cs_eS = cs_ems(indexS,2);

N0 = 1.1e25; %m^-3
f = 1/Area;
gamma = 1;
loss = 0.00;
IsatPa = h*freqP/cs_aP/tau21;
IsatSa = h*freqS/cs_aS/tau21;
IsatP = h*freqP/(cs_aP+cs_eP)/tau21;
IsatS = h*freqS/(cs_aS+cs_eS)/tau21;

z = 0:dz:len;                                         
Pp = zeros(1,length(z)); 
Psforward = zeros(1,length(z));
Psbackward = zeros(1,length(z));
Pp(1) = Pp0; %W                                      % initial pump power
Psforward(1) = Psforward0; %W                        % initial forward signal power
Psbackward(1) = Psbackward0; %W                      % initial backward signal power (GUESS)
%N2(z)
N2 = @(Pp,Ps) (Pp*f/IsatPa+Ps*f/IsatSa)./(1+Pp*f/IsatP+Ps*f/IsatS)*N0;
%signal gain coefficients
gammaS = @(Pp,Ps) N2(Pp,Ps)*(cs_eS+cs_aS)-cs_aS*N0;
gammaSE = @(Pp,Ps) cs_eS*N2(Pp,Ps);
%pumps absorption coefficients
gammaP = @(Pp,Ps) -1*N2(Pp,Ps)*(cs_aP+cs_eP)+cs_aP*N0;

%pump power
dPp_dz = @(z,Pp,Psf,Psb) -1*gammaP(Pp,Psf+Psb)*Pp;  
%forward signal power
dPsforward_dz = @(z,Pp,Psf,Psb) gammaS(Pp,Psf+Psb)*Psf+gammaSE(Pp,Psf+Psb)*2*h*freqS*dvS;
%backward signal power
dPsbackward_dz = @(z,Pp,Psf,Psb) -1*gammaS(Pp,Psf+Psb)*Psb-gammaSE(Pp,Psf+Psb)*2*h*freqS*dvS;

iterations = 0;

while iterations < maxIterations
    
    %forward propagation
    dz = abs(dz);
    for i=1:(length(z)-1)                              
        k_1 = dPp_dz(z(i),Pp(i),Psforward(i),Psbackward(i));
        j_1 = dPsforward_dz(z(i),Pp(i),Psforward(i),Psbackward(i));
        l_1 = dPsbackward_dz(z(i),Pp(i),Psforward(i),Psbackward(i));

        k_2 = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*k_1,Psforward(i)+0.5*dz*j_1,Psbackward(i)+0.5*dz*l_1);
        j_2 = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*k_1,Psforward(i)+0.5*dz*j_1,Psbackward(i)+0.5*dz*l_1);
        l_2 = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*k_1,Psforward(i)+0.5*dz*j_1,Psbackward(i)+0.5*dz*l_1);

        k_3 = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*k_2,Psforward(i)+0.5*dz*j_2,Psbackward(i)+0.5*dz*l_2);
        j_3 = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*k_2,Psforward(i)+0.5*dz*j_2,Psbackward(i)+0.5*dz*l_2);
        l_3 = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*k_2,Psforward(i)+0.5*dz*j_2,Psbackward(i)+0.5*dz*l_2);

        k_4 = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*k_3,Psforward(i)+0.5*dz*j_3,Psbackward(i)+0.5*dz*l_3);
        j_4 = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*k_3,Psforward(i)+0.5*dz*j_3,Psbackward(i)+0.5*dz*l_3);
        l_4 = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*k_3,Psforward(i)+0.5*dz*j_3,Psbackward(i)+0.5*dz*l_3);

        Pp(i+1) = Pp(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dz;
        Psforward(i+1) = Psforward(i) + (1/6)*(j_1+2*j_2+2*j_3+j_4)*dz;  % main equation
        Psbackward(i+1) = Psbackward(i) + (1/6)*(l_1+2*l_2+2*l_3+l_4)*dz;  % main equation
    end

    %check/modify boundary conditions
    if abs(Psbackward(length(z))) < maxError
        break
    else
        Psbackward(length(z)) = 0;
    end

    %backward propagation
    dz = -abs(dz);
    for i=flip(2:(length(z)))                              
        k_1 = dPp_dz(z(i),Pp(i),Psforward(i),Psbackward(i));
        j_1 = dPsforward_dz(z(i),Pp(i),Psforward(i),Psbackward(i));
        l_1 = dPsbackward_dz(z(i),Pp(i),Psforward(i),Psbackward(i));

        k_2 = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*k_1,Psforward(i)+0.5*dz*j_1,Psbackward(i)+0.5*dz*l_1);
        j_2 = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*k_1,Psforward(i)+0.5*dz*j_1,Psbackward(i)+0.5*dz*l_1);
        l_2 = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*k_1,Psforward(i)+0.5*dz*j_1,Psbackward(i)+0.5*dz*l_1);

        k_3 = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*k_2,Psforward(i)+0.5*dz*j_2,Psbackward(i)+0.5*dz*l_2);
        j_3 = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*k_2,Psforward(i)+0.5*dz*j_2,Psbackward(i)+0.5*dz*l_2);
        l_3 = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*k_2,Psforward(i)+0.5*dz*j_2,Psbackward(i)+0.5*dz*l_2);

        k_4 = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*k_3,Psforward(i)+0.5*dz*j_3,Psbackward(i)+0.5*dz*l_3);
        j_4 = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*k_3,Psforward(i)+0.5*dz*j_3,Psbackward(i)+0.5*dz*l_3);
        l_4 = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*k_3,Psforward(i)+0.5*dz*j_3,Psbackward(i)+0.5*dz*l_3);

        Pp(i-1) = Pp(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dz;
        Psforward(i-1) = Psforward(i) + (1/6)*(j_1+2*j_2+2*j_3+j_4)*dz;  % main equation
        Psbackward(i-1) = Psbackward(i) + (1/6)*(l_1+2*l_2+2*l_3+l_4)*dz;  % main equation
    end
    
    %check/modify boundary conditions
    if abs(Psforward(1)) < maxError && abs(Pp(1)-Pp0) < maxError
        break
    else
        Psforward(1) = 0;
        Pp(1) = Pp0;
    end
    
    iterations = iterations+1;
end


%calculate gain
gain = gammaS(Pp,Psforward+Psbackward);
%normalize gain
gainNorm = gain/max(gain)*Pp0;

%graphs
figure(1)
plot(z,Pp);
xlabel('z (m)');
ylabel('Pp(z) (W)');
title('Change in pump Power along fiber');

figure(2)
plot(z,Psforward);
xlabel('z (m)');
ylabel('Psforward(z) (W)');
title('Change in forward signal along fiber');

figure(3)
plot(z,Psbackward);
xlabel('z (m)');
ylabel('Psbackward(z) (W)');
title('Change in backward signal along fiber');

figure(4)
grid on
plot(z,gain);
xlabel('z (m)')
ylabel('gain (1/m)');
title('Gain along fiber');

figure(5)
grid on
plot(z,Pp,z,Psforward,z,Psbackward,z,gainNorm);
xlabel('z (m)');
ylabel('Power (W)');
title('Change in Power along fiber');
legend('pump','forward signal','backward signal','Normalized gain');

% %ODEs = @(z,Pp,Psf,Psb) [dPp_dz(z,Pp,Psf,Psb);dPsforward_dz(z,Pp,Psf,Psb);dPsbackward_dz(z,Pp,Psf,Psb)];
% ODEs = @(z,P) [dPp_dz(z,P(1),P(2),P(3));dPsforward_dz(z,P(1),P(2),P(3));dPsbackward_dz(z,P(1),P(2),P(3))];
% 
% [z, sols] = ode45(ODEs,[0 len],[Pp0 Psforward0 Psbackward0]);
% 
% %calculate gain
% gain = gammaS(sols(:,1),sols(:,2)+sols(:,3));
% %normalize gain
% gainNorm = gain/max(gain)*Pp0;
% 
% %graphs
% figure(4)
% grid on
% plot(z,sols(:,1));
% xlabel('z (m)');
% ylabel('Pp(z) (W)');
% title('Change in pump Power along fiber');
% 
% figure(5)
% grid on
% plot(z,sols(:,2));
% xlabel('z (m)');
% ylabel('Psforward(z) (W)');
% title('Change in forward signal along fiber');
% 
% figure(6)
% grid on
% plot(z,sols(:,3));
% xlabel('z (m)');
% ylabel('Psbackward(z) (W)');
% title('Change in backward signal along fiber');
% 
% figure(7)
% grid on
% plot(z,gain);
% xlabel('z (m)')
% ylabel('gain (1/m)');
% title('Gain along fiber');
% 
% figure(8)
% grid on
% plot(z,sols(:,1),z,sols(:,2),z,sols(:,3),z,gainNorm);
% xlabel('z (m)');
% ylabel('Power (W)');
% title('Change in Power along fiber');
% legend('pump','forward signal','backward signal','Normalized gain');



