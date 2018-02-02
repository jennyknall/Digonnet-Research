%numerical result for change in power across fiber for one pump and
%broadboand signal. Automatically converges. Does NOT approximate intensity in
%fiber as P(z)/Area.
%Automatically discards signal wavelengths for which gain is too low
%includes non-radiative term and loss

%close all;
clear all;

len = 3;
dz = len/300;
Pp0 = 300;
maxIterations = 401;

%ASE signal error
maxErrorS = 1e-8;

%laser sigal error
maxErrorL = 1e-4; 
errorStep = 0.005;

%cavity parameters
R1 = 1;
R2 = .01;
t1 = 1-R1;
t2 = 1-R2;

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

core_radius = 15e-6; 
Area = pi*core_radius^2; %m^2, pump area
tauRad = 0.85e-3;
tauNonRad = 1e8;
tauRatio = tauRad/tauNonRad;
tau = tauRad*tauNonRad/(tauRad+tauNonRad);
c = 3e8;
h = 6.63e-34; %J*s
kT = 4.11e-21; %J
sb = 5.67e-8; %W/m^2/K^4

%cross sectional areas 
wavelengths = (850:1150)*1e-9;
%wavelengths = (850:1600)*1e-9;
cs_absRAW = xlsread('LAS_Yb_06_02_abs.xlsx');
cs_abs = [wavelengths; interp1(cs_absRAW(:,1)*1e-9,cs_absRAW(:,2),wavelengths)].';
cs_abs(isnan(cs_abs)) = 0;
cs_emsRAW = xlsread('LAS_Yb_06_02_emi.xlsx');
cs_ems = [wavelengths; interp1(cs_emsRAW(:,1)*1e-9,cs_emsRAW(:,2),wavelengths)].';
cs_ems(isnan(cs_ems)) = 0;

%calculate mean flourecent wavelength using Mina's formula
integral1 = cs_emsRAW(:,2)./cs_emsRAW(:,1).^4;%cs_emsRAW(:,1).*cs_emsRAW(:,2)./cs_emsRAW(:,1).^5;
integral2 = cs_emsRAW(:,2)./cs_emsRAW(:,1).^5;
lambdaF_SE = trapz(cs_emsRAW(:,1),integral1)/trapz(cs_emsRAW(:,1),integral2)*1e-9;
freqF_SE = c/lambdaF_SE;

%create a variable to help me compare cross section values
crossSection = cs_abs;
crossSection(:,3) = cs_ems(:,2);
crossSection(:,4) = cs_abs(:,2)./cs_ems(:,2);

%waveleghth of laser/cooling pump
dlam = 5e-9; %resolution of the signal spectrum
lambdaP = 1030e-9;
lambdaL = 1050e-9;
lambdaS = [1035e-9:dlam:lambdaL-dlam lambdaL lambdaL+dlam:dlam:1150e-9]; %%MAYBE CHANGE THIS -- to make spacing even
freqP = c/lambdaP;
freqS = c./lambdaS;
dvS = c./(lambdaS-dlam/2)-c./(lambdaS+dlam/2); %zeros(1,length(lambdaS));%
%cross sectional areas for a given wavelength
indexP = find(round(cs_abs(:,1)*1e9 - lambdaP*1e9) == 0);
indexS = NaN(size(lambdaS));
indexL = find(lambdaS == lambdaL);
for i = 1:length(lambdaS)
    indexS(i) = find(round(cs_abs(:,1)*1e9 - lambdaS(i)*1e9) == 0);
end
cs_aP = cs_abs(indexP,2);
cs_eP = cs_ems(indexP,2);
cs_aS = cs_abs(indexS,2).';
cs_eS = cs_ems(indexS,2).';


lam = 1000e-9;
NA = .13; 
%V = 2*pi*core_radius/lam*NA;
V = 2.44;
w = core_radius*(0.65+1.619/V^1.5+2.879/V^6);
etta = 1-exp(-2*core_radius^2/w^2);

N0_wtPercent = 2.71*5;%[0.00005:.00005:.001 .001:.005:.1 .1:.1:3 3:15];
N0 = N0_wtPercent/2.71*1.81e26;  %m^-3
f = 1/Area;
gamma = 1;
totLoss = 0.0; %dB/m
loss_b = totLoss/etta*log(10)/10; %m^-1
loss_bs = loss_b/4;%m^-1
loss_ba = 3*loss_b/4;%0; %m^-1
IsatPa = h*freqP/cs_aP/tauRad;
PsatPa = IsatPa*pi*w^2/2;
IsatSa = h*freqS./cs_aS/tauRad;
IsatP = h*freqP/(cs_aP+cs_eP)/tauRad;
PsatP = IsatP*pi*w^2/2;
IsatS = h*freqS./(cs_aS+cs_eS)/tauRad;
PsatS = IsatS*pi*w^2/2;
Am = pi*w^2/2;

maxLoss = etta*cs_aP*N0*(tau/tauRad*freqF_SE/freqP-1)*10/log(10); %dB/m

% V = 1:0.01:10;
% w = core_radius*(0.65+1.619./V.^1.5+2.879./V.^6);
% PsatP = IsatP*pi*w.^2/2;
% Nup = 1-exp(-2*core_radius^2./w.^2);

%Nup = 0.01:0.0001:0.9999;
% 
%Nup = 0.9999999999
%Nup = 0.9997;
% loss_ba = .0000000000000000001;

QmaxNoLoss = -h*(tauRad/tau*freqP-freqF_SE)/(cs_aP+cs_eP)*Am*cs_aP*N0/tauRad*(-2*core_radius^2/w^2);

Popt = IsatP*pi/2*(-2*core_radius^2./log(1-etta)).*(-1*(2-etta)+sqrt(etta.^2+4*(1-etta)*cs_aP*N0/loss_ba*(tau/tauRad*lambdaP/lambdaF_SE-1)))/2./(1-etta);

%Pp0 = Popt;

%Gaussian mode intensity function
fr = @(r) 2/pi/w^2*exp(-2*r.^2/w^2);
%N2(z) - Pstot is 1xlength(lambdaS) array where each element is the sum of
%the backward and forward power for that wavelength at position 
N2 = @(Pp,Pstot,r) (Pp*fr(r)/IsatPa+sum(Pstot*fr(r)./IsatSa))/(tauRatio+1+Pp*fr(r)/IsatP+sum(Pstot*fr(r)./IsatS))*N0; 
%signal gain coefficients
gammaS = @(Pp,Pstot,r) N2(Pp,Pstot,r)*(cs_eS+cs_aS)-cs_aS*N0; %returns and array of size 1xlength(lambdaS) 
gammaSE = @(Pp,Pstot,r) cs_eS*N2(Pp,Pstot,r); %returns and array of size 1xlength(lambdaS)

comTerm1 = @(pump,Psf,Psb) (pump/PsatPa+sum((Psf+Psb)./IsatSa)/Am);
comTerm2 = @(pump,Psf,Psb) (pump/PsatP+sum((Psf+Psb)./IsatS)/Am);
%pump power - Psf and Psb are 1xlengthlambdaS arrays where each element
%is the forward and backward power (respectively) for that wavelegth and
%postion.
dPp_dz = @(z,pump,Psf,Psb) pump*cs_aP*N0/comTerm2(pump,Psf,Psb)...
    *(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb)))...
    +pump*N0*(cs_aP*sum((Psf+Psb)./IsatS)/Am-(cs_aP+cs_eP)*sum((Psf+Psb)./IsatSa)/Am)/comTerm2(pump,Psf,Psb)...
    *((exp(-2*core_radius^2/w^2)-1)-(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb)))/comTerm2(pump,Psf,Psb))...
    -loss_b*pump*(1-exp(-2*core_radius^2/w^2));
%forward signal power
dPsforward_dz = @(z,pump,Psf,Psb) -Psf.*(cs_eS+cs_aS)*N0*comTerm1(pump,Psf,Psb)/comTerm2(pump,Psf,Psb)*((exp(-2*core_radius^2/w^2)-1)...
    -(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb)))/comTerm2(pump,Psf,Psb))...
    -Psf.*cs_aS*N0*(1-exp(-2*core_radius^2/w^2))-2*h*freqS.*dvS.*cs_eS*N0*comTerm1(pump,Psf,Psb)/comTerm2(pump,Psf,Psb)...
    *((exp(-2*core_radius^2/w^2)-1)-(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb)))/comTerm2(pump,Psf,Psb))...
    -loss_b*Psf*(1-exp(-2*core_radius^2/w^2));
%backward signal power
dPsbackward_dz = @(z,pump,Psf,Psb) Psb.*(cs_eS+cs_aS)*N0*comTerm1(pump,Psf,Psb)/comTerm2(pump,Psf,Psb)*((exp(-2*core_radius^2/w^2)-1)...
    -(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb)))/comTerm2(pump,Psf,Psb))...
    +Psb.*cs_aS*N0*(1-exp(-2*core_radius^2/w^2))+2*h*freqS.*dvS.*cs_eS*N0*comTerm1(pump,Psf,Psb)/comTerm2(pump,Psf,Psb)...
    *((exp(-2*core_radius^2/w^2)-1)-(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb)))/comTerm2(pump,Psf,Psb))...
    +loss_b*Psb*(1-exp(-2*core_radius^2/w^2));



z = 0:dz:len;  
Pp = zeros(1,length(z)); 
Psforward = zeros(length(lambdaS),length(z)); 
Psbackward = zeros(length(lambdaS),length(z)); 

Ptotsignal = zeros(1,length(Pp0)); %stores total signal power emitted from fiber for each Pp0
Pabs = zeros(1,length(Pp0)); %stores values for pump power absorbed for each Pp0
dQ_dt1 = zeros(1,length(Pp0)); %stores values for dQ/dt for each Pp0
dQz_dt2 = zeros(length(Pp0),length(z)); %stores dQ/dt as a function of z for each Pp0
dQz_dtSE = zeros(length(Pp0),length(z)); %stores the SE contribution to dQ/dt as a function of z for each Pp0
dQzcorr_dtSE = zeros(length(Pp0),length(z)); %stores the correction term the SE cooling as a funciton of z for each Pp0
dQz_dtPhononLoss = zeros(length(Pp0),length(z)); %stores the phonon loss contribution to dQ/dt as a function of z for each Pp0
dQz_dtScat = zeros(length(Pp0),length(z)); %stores the scattering loss contribution to dQ/dt as a function of z for each Pp0
dSF_dz = zeros(length(Pp0),length(z)); %stores the forward signal contribution to dQ/dt as a function of z for each Pp0
dSB_dz = zeros(length(Pp0),length(z)); %stores the backward signal contribution to dQ/dt as a function of z for each Pp0 
dP_dz = zeros(length(Pp0),length(z));%stores the Pump abs contribution to dQ/dt as a function of z for each Pp0
dQ_dt2 = zeros(1,length(Pp0)); %stores values for dQ/dt for each Pp0
dQ_dtSE = zeros(1,length(Pp0)); %stores the SE contribution to dQ/dt for each Pp0
dQ_dtSEcorr = zeros(1,length(Pp0)); %stores the SE correction term to dQ/dt for each Pp0
dQ_dtPhononLoss = zeros(1,length(Pp0)); %stores the Phonon loss contribution to dQ/dt for each Pp0
dQ_dtScat = zeros(1,length(Pp0)); %stores the scattering loss contribution to dQ/dt for each Pp0
dQ_dtSF = zeros(1,length(Pp0)); %stores the forward signal contribution to dQ/dt for each Pp0
dQ_dtSB = zeros(1,length(Pp0)); %stores the backward signal contribution to dQ/dt for each Pp0
dQ_dtP = zeros(1,length(Pp0)); %stores the pump abs contribution to dQ/dt for each Pp0
maxdQz_dt2 = zeros(2,length(Pp0)); %stores the coordinates of the maximum value of dQz_dt2
                                    %first row: max values
                                    %second row: index of z where max value is
lambdaF_SL = zeros(1,length(Pp0)); %stores the mean flouresnce wavelength for each pump power calculated from the 
                                    %populations of each sublevel.
dTz = zeros(length(Pp0),length(z));  %stores the change in temp at a location z along the fiber
dTzVac = zeros(length(Pp0),length(z));  %stores the change in temp at a location z along the fiber in a vaccuum
dTavg = zeros(1,length(Pp0)); %stores the average change in temperature across the fiber 
dTVacAvg = zeros(1,length(Pp0)); %stores the average change in temperature across the fiber in a vaccuum. 
dTVacMax = zeros(1,length(Pp0)); %stores the maximum temperature change across the fiber in a vac

for p = 1:length(Pp0)
    
errorf = NaN(1,maxIterations);
errorb = NaN(1,maxIterations);
forward = NaN(1,maxIterations);
    
%calculate initial guess for Psforward(1)
nps = Area/gamma/tauRad/(cs_aP+cs_eP);
nss = Area/gamma/tauRad/(cs_aS(indexL)+cs_eS(indexL));
np0 = Pp0(p)/h/freqP;
Psforward(indexL,1) = ((gamma*cs_aS(indexL)*N0*len-log(R1*R2)/2)*nss-np0*(1-exp(-gamma*cs_aP*N0*len)*exp((gamma*cs_aS(indexL)*N0*len-log(R1*R2)/2)*nss/nps)))...
    /(1-1/sqrt(R1*R2)+(sqrt(R1*R2)-1)/R1)*h*freqS(indexL);

    
Psbackward(indexL,1) = Psforward(indexL,1)/R1;
    
allWavelengths = find(ones(1,length(lambdaS)));
wavelengthsASE = find(allWavelengths ~= indexL);
    
Pp(1) = Pp0(p); %W  

iterations = 1;
tries = 0;


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
        if sum(Psforward(:,i+1)<zeros(length(lambdaS),1)) > 0
            Psforward(Psforward(:,i+1)<zeros(length(lambdaS),1),i+1) = 0;
        end
        Psbackward(:,i+1) = Psbackward(:,i) + (1/6)*(RK1b+2*RK2b+2*RK3b+RK4b).'*dz;  % main equation
        if sum(Psbackward(:,i+1)<zeros(length(lambdaS),1)) > 0
            Psbackward(Psbackward(:,i+1)<zeros(length(lambdaS),1),i+1) = 0;
        end
    end
    
    iterations = iterations+1;

    %check/modify boundary conditions
    errorf(iterations-1) = Psforward(indexL,end)*R2/Psbackward(indexL,end)-1;
    if (sum(abs(Psbackward(wavelengthsASE,length(z))) > maxErrorS) == 0 && (sum(Psbackward(wavelengthsASE,length(z))) > 0 || sum(dvS)==0))...
            && abs(errorf(iterations-1)) < maxErrorL || iterations > maxIterations
        break
    else
        forward(iterations-1) = Psforward(indexL,1);
        Psbackward(wavelengthsASE,length(z)) = 0;
        %Psbackward(indexL,end) = Psforward(indexL,end)*R2;
        if Pp(end) > Pp(1)
            Pp(end) = Pp0(p)*exp(-cs_aP*N0*len*(1-exp(-2*core_radius^2/w^2)));
        end
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
        if sum(Psforward(:,i-1)<zeros(length(lambdaS),1)) > 0
            Psforward(Psforward(:,i-1)<zeros(length(lambdaS),1),i-1) = 0;
        end
        Psbackward(:,i-1) = Psbackward(:,i) + (1/6)*(RK1b+2*RK2b+2*RK3b+RK4b).'*dz;  % main equation
        if sum(Psbackward(:,i-1)<zeros(length(lambdaS),1)) > 0
            Psbackward(Psbackward(:,i-1)<zeros(length(lambdaS),1),i-1) = 0;
        end
    end
    
    iterations = iterations+1;
    
    %check/modify boundary conditions -- ALSO CHECK LASER CONDITIONS!!
    errorb(iterations-1) = Psbackward(indexL,1)*R1/Psforward(indexL,end)-1;
    if (sum(abs(Psforward(1)) > maxErrorS) == 0 && abs(Pp(1)-Pp0(p)) < maxErrorS)...
            && abs(errorb(iterations-1)) < maxErrorL || iterations > maxIterations
        break
    end
    
    Psforward(wavelengthsASE,1) = 0;
    Pp(1) = Pp0(p);
     
    %forwardLaser = Psforward(indexL,1);
    %Psforward(indexL,1) = forwardLaser+forwardLaser*errorStep*errorf(iterations-2)/abs(errorf(iterations-2)); 
    Psforward(indexL,1) = forward(iterations-2)+forward(iterations-2)*errorStep*errorf(iterations-2)/abs(errorf(iterations-2)); 
    Psbackward(indexL,1) = Psforward(indexL,1)/R1;
    
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
lambdaFF_ASE = 0;
lambdaFB_ASE = 0;
for i = 1:numWaves
    lambdaFF_ASE = lambdaFF_ASE+lambdaS(i)*Psforward(i,length(z));
    lambdaFB_ASE = lambdaFB_ASE+lambdaS(i)*Psbackward(i,1);
end
lambdaFF_ASE = lambdaFF_ASE/sum(Psforward(:,length(z)));
lambdaFB_ASE = lambdaFB_ASE/sum(Psbackward(:,1));


%calculate dQ/dt
%method 2: solving differential equation
T1 = @(Ps) sum(Ps./IsatSa)*h*freqF_SE/tauRad;%sum(Ps.*cs_aS);%
T2 = @(Ps) sum(Ps./IsatS)*h*freqF_SE/tauRad;%sum(Ps.*(cs_aS+cs_eS));%
E = sum(cs_eS.*freqS.*dvS);
%find dQ/dt as a function of z
for zz = 1:length(z)
    dP_dz(p,zz) = -N0*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')...
        *Pp(zz)*(cs_aP+cs_eP)*((exp(-2*core_radius^2/w^2)-1)...
        -(1+tauRatio)*log((tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')))...
        /comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).'))...
        -N0*Pp(zz)*cs_aP*(1-exp(-2*core_radius^2/w^2));
    
    dSF_dz(p,zz) = -N0*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')...
        *(T2(Psforward(:,zz).')+2*h*sum(cs_eS.*freqS.*dvS))*((exp(-2*core_radius^2/w^2)-1)...
        -(1+tauRatio)*log((tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')))...
        /comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).'))...
        -N0*T1(Psforward(:,zz).')*(1-exp(-2*core_radius^2/w^2));

    dSB_dz(p,zz) = N0*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')...
        *(T2(Psbackward(:,zz).')+2*h*sum(cs_eS.*freqS.*dvS))*((exp(-2*core_radius^2/w^2)-1)...
        -(1+tauRatio)*log((tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')))...
        /comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).'))...
        +N0*T1(Psbackward(:,zz).')*(1-exp(-2*core_radius^2/w^2));

    dQz_dtSE(p,zz) = -Am*N0*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')*h*freqF_SE/tauRad...
        *log((tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')));

    dQzcorr_dtSE(p,zz) = -4*h*N0*E*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')...
        *((exp(-2*core_radius^2/w^2)-1)...
        -(1+tauRatio)*log((tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).')))...
        /comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).'));
    
    dQz_dtPhononLoss(p,zz) = (Pp(zz)+sum(Psforward(:,zz)+Psbackward(:,zz)))*loss_ba*(1-exp(-2*core_radius^2/w^2));
    
    dQz_dtScat(p,zz) = (Pp(zz)+sum(Psforward(:,zz)+Psbackward(:,zz)))*loss_bs*(1-exp(-2*core_radius^2/w^2));
    
    dQz_dt2(p,zz) = -dP_dz(p,zz)-dSF_dz(p,zz)+dSB_dz(p,zz)-dQz_dtSE(p,zz)+dQzcorr_dtSE(p,zz)+dQz_dtPhononLoss(p,zz);
end
[maxdQz_dt2(1,p), maxdQz_dt2(2,p)] = min(dQz_dt2(p,:));
%integrate over z to find dQ/dt for whole fiber
dQ_dt2(p) = trapz(z,dQz_dt2(p,:));
dQ_dtSE(p) = trapz(z,-dQz_dtSE(p,:));
dQ_dtSEcorr(p) = trapz(z,dQzcorr_dtSE(p,:));
dQ_dtPhononLoss(p) = trapz(z,dQz_dtPhononLoss(p,:));
dQ_dtScat(p) = trapz(z,-dQz_dtScat(p,:));
dQ_dtSF(p) = trapz(z,-dSF_dz(p,:));
dQ_dtSB(p) = trapz(z,dSB_dz(p,:));
dQ_dtP(p) = trapz(z,-dP_dz(p,:));
%method 1: Pabs-Ps 
Ptotsignal(p) = sum(Psforward(:,length(z))-Psforward(:,1)+Psbackward(:,1)-Psbackward(:,length(z)));
Pabs(p) = Pp(1)-Pp(length(z));
dQ_dt1(p) = Pabs(p)-Ptotsignal(p)+dQ_dtSE(p)+dQ_dtSEcorr(p)+dQ_dtScat(p);

%calculate change in temperature
b = 67.3e-6; %m
%b = 125e-6; %m
SA = 2*pi*b*len;
density = 2.65e3; %kg/m^3
cv = 680; %J/kg/K
mcv = pi*b^2*len*density*cv;
rTemp = 300;
emis_f = 1;
emis_c = 0.05;
radius_c = 2.5e-3; %m
chi = (1-emis_c)*emis_f*b/emis_c/radius_c;
thermalC = 81.4; %W/m^2/K
dTz(p,:) = dQz_dt2(p,:)/2/pi/b/thermalC; %K
dTzVac(p,:) = (rTemp^4+dQz_dt2(p,:)*(1+chi)/(sb*2*pi*b*emis_f)).^0.25;
dTVacAvg(p) = mean(dTzVac(p,:));
dTVacMax(p) = min(dTzVac(p,:));
dTavg(p) = mean(dTz(p,:)); %K

TmaxNoLoss = (rTemp^4+QmaxNoLoss*(1+chi)/(sb*2*pi*b*emis_f)).^0.25;


%calculate N2(z)
N2z = zeros(1,length(z));
N1z = zeros(1,length(z));
for r = 0:dr:core_radius-dr
    for i = 1:length(z)
        N2z(i) = N2z(i)+N2(Pp(i),Psforward(:,i).'+Psbackward(:,i).',r)/numPoints;
    end
end
N1z = N0*ones(1,length(z))-N2z;

%Calculating N2j and N1i 
N2j = zeros(g2,length(z)); %dim-> height:each sublevel
                            %      length:population of each sublevel at each position

%find N2/N21 and N1/N11
N2_N21 = 0; 
for j = 1:g2
    delE = (E2(j)-E21)*h*c;
    N2_N21 = N2_N21 + exp(-delE/kT);
end
N1_N11 = 0;
for i = 1:g1
    delE = (E1(i)-E11)*h*c;
    N1_N11 = N1_N11 + exp(-delE/kT);
end

for zz = 1:length(z) %at each position, find the populations of each sublevel.
    for j = 1:g2
        delE = (E2(j)-E21)*h*c;
        N2j(j,zz) = exp(-delE/kT)*N2z(zz)/N2_N21;
    end
end

%calculate mean flourecence wavelength along z using sublevel populations
lambdaFz_SL = zeros(1,length(z));
for zz = 1:length(z)
    sum1 =g1*N2j(1,zz)+g1*N2j(2,zz)+g1*N2j(3,zz);
    sum2 = N2j(1,zz)*(E21-E11)+N2j(1,zz)*(E21-E12)+N2j(1,zz)*(E21-E13)+N2j(1,zz)*(E21-E14)+...
        N2j(2,zz)*(E22-E11)+N2j(2,zz)*(E22-E12)+N2j(2,zz)*(E22-E13)+N2j(2,zz)*(E22-E14)+...
        N2j(3,zz)*(E23-E11)+N2j(3,zz)*(E23-E12)+N2j(3,zz)*(E23-E13)+N2j(3,zz)*(E23-E14);
    lambdaFz_SL(zz) = sum1/sum2;
end
lambdaF_SL(p) = mean(lambdaFz_SL);


% %graphs
%Pump power
figure(1)
grid on
hold on
plot(z,Pp);
xlabel('z (m)');
ylabel('Pp(z) (W)');
title('Change in pump Power along fiber');

%forward signal
figure(2)
hold on
grid on
for i = wavelengthsASE
    plot(z,Psforward(i,:),'DisplayName',sprintf('%g nm',lambdaS(i)*1e9));
end
xlabel('z (m)');
ylabel('Forward ASE Power (W)');
title('Change in Forward ASE Signal Along The Fiber');
legend('show');

%backward signal
figure(3)
hold on
grid on
for i = wavelengthsASE
    plot(z,Psbackward(i,:),'DisplayName',sprintf('%g nm',lambdaS(i)*1e9));
end
xlabel('z (m)');
ylabel('Backward ASE Power (W)');
title('Change in Backward ASE Signal Along The Fiber');

%Total ASE power
figure(4)
hold on 
grid on
if length(lambdaS) > 1
    plot(z,sum(Psforward(wavelengthsASE,:)));
    plot(z,sum(Psbackward(wavelengthsASE,:)));
else
    plot(z,Psforward);
    plot(z,Psbackward);
end
xlabel('z (m)');
ylabel('Power (W)');
title('ASE Power along fiber');
legend('forward signal','backward signal');

%pump and laser signal
figure(5)
hold on 
grid on
plot(z,Pp);
plot(z,Psforward(indexL,:));
plot(z,Psbackward(indexL,:));
xlabel('z (m)');
ylabel('Power (W)');
title('Laser Power along fiber');
legend('Pump','forward signal','backward signal');

%forward Signal Spectrum
figure(6)
hold on
grid on
PsforwardNorm = Psforward(wavelengthsASE,length(z))*1e3/dlam/1e9;
plot(lambdaS(wavelengthsASE)*1e9,PsforwardNorm,'DisplayName',sprintf('%g mW, mean = %g nm',Pp0(p)*1e3,lambdaFF_ASE*1e9));
xlabel('wavelength (nm)');
ylabel('power (mW/nm)');
title('Normalized Spectrum of forward signal')

%backward Signal Spectrum
figure(7)
hold on
grid on
PsbackwardNorm = Psbackward(wavelengthsASE,1)*1e3/dlam/1e9;
plot(lambdaS(wavelengthsASE)*1e9,PsbackwardNorm,'DisplayName',sprintf('%g mW, mean = %g nm',Pp0(p)*1e3,lambdaFB_ASE*1e9));
xlabel('wavelength (nm)');
ylabel('power (mW/nm)');
title('Normalized Spectrum of backward signal')

% %Signal Gain
% figure(8) 
% hold on
% grid on
% for i = 95%1:length(lambdaS)
%     plot(z,gain(i,:),'DisplayName',sprintf('%g nm',lambdaS(i)*1e9));
% end
% xlabel('z (m)');
% ylabel('gain(z) (1/m)');
% title('Signal gain along fiber');
% legend('show')
% 
% %Signal Gain vs. Wavelength
% figure(9)
% hold on
% grid on
% plot(lambdaS*1e9,avgGain)
% plot(lambdaS*1e9,maxGain)
% xlabel('wavelength (nm)');
% ylabel('average gain(z) (1/m)');
% title('Average/Max signal gain along fiber vs. wavelength');
% legend('average gain','max gain');

%N2(z)
figure(10)
grid on
hold on
plot(z,N2z/N0,'DisplayName',sprintf('%g mW',Pp0(p)*1e3));
xlabel('z (m)');
ylabel('N2(z)/N0');
title('N2 along fiber');

%dQ/dt as a function of z
figure(11)
%grid on 
box on
hold on
plot(z/len,dQz_dt2(p,:)*1e3,'DisplayName',sprintf('%g mW',Pp0(p)*1e3));
%xlabel('Position z along the fiber (m)');
xlabel('Normalized position z along the fiber (%)');
ylabel('Heat extracted per unit length (mW/m)');
title('dQ/dt along the fiber');

%convergence error
figure(12)
hold on
plot(1:sum(isfinite(errorf)),errorf(isfinite(errorf)));
plot(1:sum(isfinite(errorb)),errorb(isfinite(errorb)));
xlabel('iteration');
ylabel('error')
title('error vs. iteration')
grid on
legend('forward error','backward error');

%change in forward signal with convergence
figure(13)
hold on
plot(1:sum(isfinite(forward)),forward(isfinite(forward)));
xlabel('iteration');
ylabel('forward signal')
title('forward signal vs. iteration')
grid on

% % 
% %dT as a function of z
% figure(14)
% hold on
% grid on
% plot(z,dTz(p,:),'DisplayName',sprintf('P_{pump} = %g mW',Pp0(p)*1e3));
% xlabel('z (m)');
% ylabel('Change In Temperature (K)');
% title('Change in Temperature Along the Fiber in Air');

end

% figure(5) %forward signal spectrum
% legend('show')
% figure(6) %backward signal spectrum
% legend('show')
% figure(8) %N2(z)
% legend('show');
% figure(13) %dTz
% legend('show');
% figure(9) %dQ/dt(z)
% legend('show')
% plot(z(maxdQz_dt2(2,:)),maxdQz_dt2(1,:));


