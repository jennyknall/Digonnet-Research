%numerical result for change in power across fiber for one pump and
%broadboand signal. Automatically converges. Does NOT approximate intensity in
%fiber as P(z)/Area

close all;
clear all;

len = 10;
dz = .01;
Pp0 = .5;%[.01:.01:1 2 3];%.01:.01:.1;%.02:.01:.2;%.001:.01:.1;%[.01:.01:1];
%Psforward0 = 0;
%Psbackward0 = ones(21,1)*Pp0/21; %%GUESS THIS!!!!
maxIterations = 10;
maxError = 1e-7;

core_radius = 1.5e-6; 
Area = pi*core_radius^2; %m^2, pump area
tauRad = 1e-3;
c = 3e8;
h = 6.63e-34; %J*s

%cross sectional areas 
wavelengths = (870:1050)*1e-9;
cs_absRAW = xlsread('abs_ZBLAN.xlsx');
cs_abs = [wavelengths; interp1(cs_absRAW(:,1),cs_absRAW(:,2),wavelengths)*1e-24].';
cs_emsRAW = xlsread('emm_ZBLAN.xlsx');
cs_ems = [wavelengths; interp1(cs_emsRAW(:,1),cs_emsRAW(:,2),wavelengths)*1e-24].';

%calculate mean flourecent wavelength using Mina's formula
integral1 = cs_emsRAW(:,1).*cs_emsRAW(:,2)./cs_emsRAW(:,1).^5;
integral2 = cs_emsRAW(:,2)./cs_emsRAW(:,1).^5;
lambda_f = trapz(cs_emsRAW(:,1),integral1)/trapz(cs_emsRAW(:,1),integral2);
freqF = c/lambda_f;

%create a variable to help me compare cross section values
crossSection = cs_abs;
crossSection(:,3) = cs_ems(:,2);
crossSection(:,4) = cs_abs(:,2)./cs_ems(:,2);

%waveleghth of laser/cooling pump
dlam = 10e-9; %resolution of the signal spectrum
lambdaP = 935e-9;
lambdaS = [930e-9:dlam:1050e-9]; 
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

ncore = 1.49;
nclad = 1.4889;
lam = 1000e-9;
NA = .205; 
V = 2*pi*core_radius/lam*NA;
w = core_radius*(0.65+1.619/V^1.5+2.879/V^6);


N0 = 1.1e25;% 1.63e26; %m^-3
f = 1/Area;
gamma = 1;
loss = 0.00;
IsatPa = h*freqP/cs_aP/tauRad;
PsatPa = IsatPa*pi*w^2/2;
IsatSa = h*freqS./cs_aS/tauRad;
IsatP = h*freqP/(cs_aP+cs_eP)/tauRad;
PsatP = IsatP*pi*w^2/2;
IsatS = h*freqS./(cs_aS+cs_eS)/tauRad;
Am = pi*w^2/2;


%Gaussian mode intensity function
fr = @(r) 2/pi/w^2*exp(-2*r.^2/w^2);
%N2(z) - Pstot is 1xlength(lambdaS) array where each element is the sum of
%the backward and forward power for that wavelength at position 
N2 = @(Pp,Pstot,r) (Pp*fr(r)/IsatPa+sum(Pstot*fr(r)./IsatSa))/(1+Pp*fr(r)/IsatP+sum(Pstot*fr(r)./IsatS))*N0; 
%signal gain coefficients
gammaS = @(Pp,Pstot,r) N2(Pp,Pstot,r)*(cs_eS+cs_aS)-cs_aS*N0; %returns and array of size 1xlength(lambdaS) 
gammaSE = @(Pp,Pstot,r) cs_eS*N2(Pp,Pstot,r); %returns and array of size 1xlength(lambdaS)

comTerm1 = @(pump,Psf,Psb) (pump/PsatPa+sum((Psf+Psb)./IsatSa)/Am);
comTerm2 = @(pump,Psf,Psb) (pump/PsatP+sum((Psf+Psb)./IsatS)/Am);
%pump power - Psf and Psb are 1xlength(lambdaS) arrays where each element
%is the forward and backward power (respectively) for that wavelegth and
%postion.
dPp_dz = @(z,pump,Psf,Psb) pump*cs_aP*N0/comTerm2(pump,Psf,Psb)...
    *log((1+comTerm2(pump,Psf,Psb)*exp(-2*core_radius^2/w^2))/(1+comTerm2(pump,Psf,Psb)))...
    +pump*N0*(cs_aP*sum((Psf+Psb)./IsatS)/Am-(cs_aP+cs_eP)*sum((Psf+Psb)./IsatSa)/Am)/comTerm2(pump,Psf,Psb)...
    *((exp(-2*core_radius^2/w^2)-1)-log((1+comTerm2(pump,Psf,Psb)*exp(-2*core_radius^2/w^2))/(1+comTerm2(pump,Psf,Psb)))/comTerm2(pump,Psf,Psb));
%forward signal power
dPsforward_dz = @(z,pump,Psf,Psb) -Psf.*(cs_eS+cs_aS)*N0*comTerm1(pump,Psf,Psb)/comTerm2(pump,Psf,Psb)*((exp(-2*core_radius^2/w^2)-1)...
    -log((1+comTerm2(pump,Psf,Psb)*exp(-2*core_radius^2/w^2))/(1+comTerm2(pump,Psf,Psb)))/comTerm2(pump,Psf,Psb))...
    -Psf.*cs_aS*N0*(1-exp(-2*core_radius^2/w^2))-2*h*freqS.*dvS.*cs_eS*N0*comTerm1(pump,Psf,Psb)/comTerm2(pump,Psf,Psb)...
    *((exp(-2*core_radius^2/w^2)-1)-log((1+comTerm2(pump,Psf,Psb)*exp(-2*core_radius^2/w^2))/(1+comTerm2(pump,Psf,Psb)))/comTerm2(pump,Psf,Psb));
%backward signal power
dPsbackward_dz = @(z,pump,Psf,Psb) Psb.*(cs_eS+cs_aS)*N0*comTerm1(pump,Psf,Psb)/comTerm2(pump,Psf,Psb)*((exp(-2*core_radius^2/w^2)-1)...
    -log((1+comTerm2(pump,Psf,Psb)*exp(-2*core_radius^2/w^2))/(1+comTerm2(pump,Psf,Psb)))/comTerm2(pump,Psf,Psb))...
    +Psb.*cs_aS*N0*(1-exp(-2*core_radius^2/w^2))+2*h*freqS.*dvS.*cs_eS*N0*comTerm1(pump,Psf,Psb)/comTerm2(pump,Psf,Psb)...
    *((exp(-2*core_radius^2/w^2)-1)-log((1+comTerm2(pump,Psf,Psb)*exp(-2*core_radius^2/w^2))/(1+comTerm2(pump,Psf,Psb)))/comTerm2(pump,Psf,Psb));


z = 0:dz:len;  
Pp = zeros(1,length(z)); 
Psforward = zeros(length(lambdaS),length(z)); 
Psbackward = zeros(length(lambdaS),length(z)); 

Ptotsignal = zeros(1,length(Pp0)); %stores total signal power emitted from fiber for each Pp0
Pabs = zeros(1,length(Pp0)); %stores values for pump power absorbed for each Pp0
dQ_dt1 = zeros(1,length(Pp0)); %stores values for dQ/dt for each Pp0
dQz_dt2 = zeros(length(Pp0),length(z)); %stores dQ/dt as a function of z for each Pp0
dQz_dtSE = zeros(length(Pp0),length(z)); %stores the SE contribution to dQ/dt as a function of z for each Pp0
dSF_dz = zeros(length(Pp0),length(z)); %stores the forward signal contribution to dQ/dt as a function of z for each Pp0
dSB_dz = zeros(length(Pp0),length(z)); %stores the backward signal contribution to dQ/dt as a function of z for each Pp0 
dP_dz = zeros(length(Pp0),length(z));%stores the Pump abs contribution to dQ/dt as a function of z for each Pp0
dQ_dt2 = zeros(1,length(Pp0)); %stores values for dQ/dt for each Pp0
dQ_dtSE = zeros(1,length(Pp0)); %stores the SE contribution to dQ/dt for each Pp0
dQ_dtSF = zeros(1,length(Pp0)); %stores the forward signal contribution to dQ/dt for each Pp0
dQ_dtSB = zeros(1,length(Pp0)); %stores the backward signal contribution to dQ/dt for each Pp0
dQ_dtP = zeros(1,length(Pp0)); %stores the pump abs contribution to dQ/dt for each Pp0
maxdQz_dt2 = zeros(2,length(Pp0)); %stores the coordinates of the maximum value of dQz_dt2
                                    %first row: max values
                                    %second row: index of z where max value is

for p = 1:length(Pp0)
    
Psbackward(:,1) = ones(length(lambdaS),1)*Pp0(p)/length(lambdaS)/2; %%GUESS THIS!!!!
Pp(1) = Pp0(p); %W  

iterations = 0;


while iterations <= maxIterations
 
    %forward propagation
    dz = abs(dz);
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
    
    %check/modify boundary conditions
    if sum(abs(Psbackward(:,length(z))) > maxError) == 0
        break
    else
        Psbackward(:,length(z)) = 0;
    end
    
    %backward propagation
    dz = -abs(dz);
    for i=flip(2:length(z)) % calculation loop

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
       
        Pp(i-1) = Pp(i) + (1/6)*(RK1p+2*RK2p+2*RK3p+RK4p)*dz;
        Psforward(:,i-1) = Psforward(:,i) + (1/6)*(RK1f+2*RK2f+2*RK3f+RK4f).'*dz;  % main equation
        Psbackward(:,i-1) = Psbackward(:,i) + (1/6)*(RK1b+2*RK2b+2*RK3b+RK4b).'*dz;  % main equation
    end
    
    iterations = iterations+1;
    
    %check/modify boundary conditions
    if (sum(abs(Psforward(:,1)) > maxError) == 0 && abs(Pp(1)-Pp0(p)) < maxError) || iterations > maxIterations
        break
    else
        Psforward(:,1) = 0;
        Pp(1) = Pp0(p);
    end
end

%calculate gain 
gain = zeros(length(lambdaS),length(z));
avgGain = zeros(1,length(lambdaS)); %stores the average gain across the fiber for each wavelength
maxGain = zeros(1,length(lambdaS)); %stores the maximum gain across the fiber for each wavelength
numPoints = 5;
dr = core_radius/numPoints;
for r = 0:dr:core_radius-dr %takes the average gain across radius of fiber
    for i = 1:length(z)
        gain(:,i) = gain(:,i)+gammaS(Pp(i),Psforward(:,i).'+Psbackward(:,i).',r).'/numPoints;
    end
end
for i = 1:length(lambdaS)
    avgGain(i) = mean(gain(i,:));
    maxGain(i) = max(gain(i,:));
end

%calculate mean flourescent wavelength of ASE spectrum
numWaves =length(lambdaS);
mlambda_f = 0;
mlambda_b = 0;
for i = 1:numWaves
    mlambda_f = mlambda_f+lambdaS(i)*Psforward(i,length(z));
    mlambda_b = mlambda_b+lambdaS(i)*Psbackward(i,1);
end
mlambda_f = mlambda_f/sum(Psforward(:,length(z)));
mlambda_b = mlambda_b/sum(Psbackward(:,1));

%calculate dQ/dt
%method 2: solving differential equation
T1 = @(Ps) sum(Ps.*cs_aS);
T2 = @(Ps) sum(Ps.*(cs_aS+cs_eS));
%find dQ/dt as a function of z
for zz = 1:length(z)
    dP_dz(p,zz) = -N0*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')...
        *Pp(zz)*(cs_aP+cs_eP)*((exp(-2*core_radius^2/w^2)-1)...
        -log((1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')*exp(-2*core_radius^2/w^2))/(1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')))...
        /comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).'))...
        -N0*Pp(zz)*cs_aP*(1-exp(-2*core_radius^2/w^2));
    
    dSF_dz(p,zz) = -N0*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')...
        *(T2(Psforward(:,zz).')+2*h*sum(cs_eS.*freqS.*dvS))*((exp(-2*core_radius^2/w^2)-1)...
        -log((1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')*exp(-2*core_radius^2/w^2))/(1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')))...
        /comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).'))...
        -N0*T1(Psforward(:,zz).')*(1-exp(-2*core_radius^2/w^2));

    dSB_dz(p,zz) = N0*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')...
        *(T2(Psbackward(:,zz).')+2*h*sum(cs_eS.*freqS.*dvS))*((exp(-2*core_radius^2/w^2)-1)...
        -log((1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')*exp(-2*core_radius^2/w^2))/(1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')))...
        /comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).'))...
        +N0*T1(Psbackward(:,zz).')*(1-exp(-2*core_radius^2/w^2));

    dQz_dtSE(p,zz) = -Am*N0*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')*h*freqF/tauRad...
        *log((1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')*exp(-2*core_radius^2/w^2))/(1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')));

     dQz_dt2(p,zz) = -dP_dz(p,zz)-dSF_dz(p,zz)+dSB_dz(p,zz)-dQz_dtSE(p,zz);
%     dQz_dt2(p,zz) = N0*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')...
%         *(Pp(zz)*(cs_aP+cs_eP)-T2(Psforward(:,zz).',Psbackward(:,zz).')-4*h*sum(cs_eS.*freqS.*dvS))...
%         *((exp(-2*core_radius^2/w^2)-1)...
%         -log((1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')*exp(-2*core_radius^2/w^2))/(1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')))...
%         /comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).'))...
%         +N0*(Pp(zz)*cs_aP-T1(Psforward(:,zz).',Psbackward(:,zz).'))*(1-exp(-2*core_radius^2/w^2))...
%         +dQz_dtSE(p,zz);
end
[maxdQz_dt2(1,p), maxdQz_dt2(2,p)] = max(dQz_dt2(p,:));
%integrate over z to find dQ/dt for whole fiber
dQ_dt2(p) = trapz(z,dQz_dt2(p,:));
dQ_dtSE(p) = trapz(z,-dQz_dtSE(p,:));
dQ_dtSF(p) = trapz(z,-dSF_dz(p,:));
dQ_dtSB(p) = trapz(z,dSB_dz(p,:));
dQ_dtP(p) = trapz(z,-dP_dz(p,:));
%method 1: Pabs-Ps 
Ptotsignal(p) = sum(Psforward(:,length(z))-Psforward(:,1)+Psbackward(:,1)-Psbackward(:,length(z)));
Pabs(p) = Pp(1)-Pp(length(z));
dQ_dt1(p) = Pabs(p)-Ptotsignal(p)+dQ_dtSE(p);


%calculate N2(z)
N2z = zeros(1,length(z));
for r = 0:dr:core_radius-dr
    for i = 1:length(z)
        N2z(i) = N2z(i)+N2(Pp(i),Psforward(:,i).'+Psbackward(:,i).',r)/numPoints;
    end
end


%graphs
%Pump power
%approx = Pp0(p)*exp(-cs_aP*N0*z*(1-exp(-2*core_radius^2/w^2)));
% figure(1)
% grid on
% hold on
% plot(z,Pp);
% %plot(z,log(approx));
% xlabel('z (m)');
% ylabel('Pp(z) (W)');
% title('Change in pump Power along fiber');
%legend('non approx','approx');

%forward signal
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

%backward signal
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

%pump and signals
figure(4)
hold on 
grid on
plot(z,Pp);
if length(lambdaS) > 1
    plot(z,sum(Psforward));
    plot(z,sum(Psbackward));
else
    plot(z,Psforward);
    plot(z,Psbackward);
end
xlabel('z (m)');
ylabel('Power (W)');
title('Change in Power along fiber');
legend('pump','forward signal','backward signal');

% %forward Signal Spectrum
% figure(5)
% hold on
% grid on
% PsforwardNorm = Psforward(:,length(z))*1e3/dlam/1e9;
% plot(lambdaS*1e9,PsforwardNorm,'DisplayName',sprintf('%g mW, mean = %g nm',Pp0(p)*1e3,mlambda_f*1e9));
% xlabel('wavelength (nm)');
% ylabel('power (mW/nm)');
% title('Normalized Spectrum of forward signal')
% 
% %backward Signal Spectrum
% figure(6)
% hold on
% grid on
% PsbackwardNorm = Psbackward(:,1)*1e3/dlam/1e9;
% plot(lambdaS*1e9,PsbackwardNorm,'DisplayName',sprintf('%g mW, mean = %g nm',Pp0(p)*1e3,mlambda_b*1e9));
% xlabel('wavelength (nm)');
% ylabel('power (mW/nm)');
% title('Normalized Spectrum of backward signal')

%Signal Gain
figure(7) 
hold on
grid on
for i = 1:length(lambdaS)
    plot(z,gain(i,:),'DisplayName',sprintf('%g nm',lambdaS(i)*1e9));
end
xlabel('z (m)');
ylabel('gain(z) (1/m)');
title('Signal gain along fiber');
legend('show')

%Signal Gain vs. Wavelength
figure(13)
hold on
grid on
plot(lambdaS*1e9,avgGain)
plot(lambdaS*1e9,maxGain)
xlabel('wavelength (nm)');
ylabel('average gain(z) (1/m)');
title('Average/Max signal gain along fiber vs. wavelength');
legend('average gain','max gain');

%N2(z)
figure(8)
grid on
hold on
plot(z,N2z);
plot(z,N0*ones(length(z)))
xlabel('z (m)');
ylabel('N2(z) (1/m^3)');
title('N2 along fiber');
% 
% %dQ/dt as a function of z
% figure(9)
% grid on 
% hold on
% plot(z,dQz_dt2(p,:),'DisplayName',sprintf('%g mW',Pp0(p)*1e3));
% xlabel('z (m)');
% ylabel('dQ/dt (J/s)');
% title('dQ/dt along the fiber');

end
% 
% figure(5) %forward signal spectrum
% legend('show')
% figure(6) %backward signal spectrum
% legend('show')
% figure(9) %dQ/dt(z)
% legend('show')
% plot(z(maxdQz_dt2(2,:)),maxdQz_dt2(1,:));
% 
% %total signal power vs. absorped pump power
% figure(10)
% grid on
% plot(Pabs*1e3,Ptotsignal*1e3)
% xlabel('Pabs (mW)');
% ylabel('total signal power emmitted (mW)');
% title('total signal power vs. absorbed pump power');

%Q
figure(11)
grid on
hold on
plot(Pp0,dQ_dt1);
plot(Pp0,dQ_dt2);
xlabel('Pp(0) (W)');
ylabel('dQ/dt (J/s)')
title('Q vs Initial pump power');
legend('Pabs-Ps','differential equation');

%SE and ASE contribution to Q
figure(12)
grid on
hold on
plot(Pp0,Ptotsignal);
plot(Pp0,dQ_dtSE);
xlabel('Pp(0) (W)');
ylabel('dQ_SE/dt (J/s)')
title('SE contribution to Q vs. Initial pump power');
legend('ASE','other SE');





% me = (Pp(length(z))-Pp(1))/len(1);
% ma = (approx(length(z))-approx(1))/len(1);
% mratio = me/ma
% logFactor = (1+exp(-2*core_radius^2/w^2))/2
% slopefactor = log(logFactor)


