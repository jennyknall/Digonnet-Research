%numerical result for change in power across fiber for one pump and
%broadboand signal. Automatically converges. Does NOT approximate intensity in
%fiber as P(z)/Area.
%Automatically discards signal wavelengths for which gain is too low
%includes non-radiative term and loss

close all;
clear all;

len = 3; %length of cavity (m)
dz = len/300;
Pp0 = 300;%40:20:60;%[.8:.2:2 1000:300:2000];%.2:.1:1;%.2:.1:.5;%.1:.1:.5;%0.06:0.02:.20;%:0.1:1;%.01:.05:1;%[.09:.03:1 2 3]; %pump power (W)
loops = 30; %humber of loops for heat distribution equilization. allowed values: 1,2,3,4,6,7,9,13

maxIterations = 400; %convergence iterations
maxError = 1e-4; 
errorStep = 0.005; 

errorf = NaN(1,maxIterations);

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

%core_radius = 2*3.1e-6; 
core_radius = 15e-6; 
Area = pi*core_radius^2; %m^2, pump area
tauRad = 0.85e-3;
tauNonRad = .176e9;
tauRatio = tauRad/tauNonRad;
c = 3e8;
h = 6.63e-34; %J*s
kT = 4.11e-21; %J

%cross sectional areas 
wavelengths = (850:1150)*1e-9;
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

%waveleghth of laser/pump
dlam = 1e-9; %linewidth of laser
lambdaP = 1030e-9; %pump wavelength
lambdaS = 1050e-9; %laser wavelength
freqP = c/lambdaP;
freqS = c/lambdaS;
dvS = 0;%c/(lambdaS-dlam/2)-c/(lambdaS+dlam/2); %zeros(1,length(lambdaS));%
%cross sectional areas for a given wavelength
indexP = find(round(cs_abs(:,1)*1e9 - lambdaP*1e9) == 0);
indexS = find(round(cs_abs(:,1)*1e9 - lambdaS*1e9) == 0);
cs_aP = cs_abs(indexP,2); 
cs_eP = cs_ems(indexP,2);
cs_aS = cs_abs(indexS,2);
cs_eS = cs_ems(indexS,2);


lam = 1000e-9;
NA = .13; 
%V = 2*pi*core_radius/lam*NA;
V = 2.44;
w = core_radius*(0.65+1.619/V^1.5+2.879/V^6);
%w = core_radius; %% use for MMF
etta = 1-exp(-2*core_radius^2/w^2);


N0_wtPercent = 2.71*5;%[0.00005:.00005:.001 .001:.005:.1 .1:.1:3 3:15];
N0 = N0_wtPercent/2.71*1.81e26;  %m^-3
%N0 = 15*1.81e26; %m^-3
Am = pi*w^2/2;
gamma = 1;
totLoss = 0.02; %dB/m
loss_b = totLoss/etta*log(10)/10; %m^-1
loss_bs = 1*loss_b/4;%m^-1
loss_ba = 3*loss_b/4;%0; %m^-1
IsatPa = h*freqP/cs_aP/tauRad;
PsatPa = IsatPa*pi*w^2/2;
IsatSa = h*freqS/cs_aS/tauRad;
IsatP = h*freqP/(cs_aP+cs_eP)/tauRad;
PsatP = IsatP*pi*w^2/2;
IsatS = h*freqS/(cs_aS+cs_eS)/tauRad;
PsatS = IsatS*Am;





%Gaussian mode intensity function
fr = @(r) 2/pi/w^2*exp(-2*r.^2/w^2);
%fr = @(r) (r<core_radius)*1/Am; %%use for MMF
%N2(z) - Pstot is the sum of the backward and forward power for that position 
N2 = @(Pp,Pstot,r) (Pp*fr(r)/IsatPa+Pstot*fr(r)/IsatSa)/(tauRatio+1+Pp*fr(r)/IsatP+Pstot*fr(r)/IsatS)*N0; 
%signal gain coefficients
gammaS = @(Pp,Pstot,r) N2(Pp,Pstot,r)*(cs_eS+cs_aS)-cs_aS*N0; %returns and array of size 1xlength(lambdaS) 
gammaSE = @(Pp,Pstot,r) cs_eS*N2(Pp,Pstot,r); %returns and array of size 1xlength(lambdaS)

comTerm1 = @(pump,Psf,Psb) (pump/PsatPa+(Psf+Psb)/IsatSa/Am);
comTerm2 = @(pump,Psf,Psb) (pump/PsatP+(Psf+Psb)/IsatS/Am);
%pump power - Psf and Psb is the value of the forward and backward power (respectively) for that postion.
dPp_dz = @(z,pump,Psf,Psb) pump*cs_aP*N0/comTerm2(pump,Psf,Psb)...
    *(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb)))...
    +pump*N0*(cs_aP*(Psf+Psb)/IsatS/Am-(cs_aP+cs_eP)*(Psf+Psb)/IsatSa/Am)/comTerm2(pump,Psf,Psb)...
    *((exp(-2*core_radius^2/w^2)-1)-(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb)))/comTerm2(pump,Psf,Psb))...
    -loss_b*pump*(1-exp(-2*core_radius^2/w^2));
%forward signal power
dPsforward_dz = @(z,pump,Psf,Psb) -Psf*(cs_eS+cs_aS)*N0*comTerm1(pump,Psf,Psb)/comTerm2(pump,Psf,Psb)*((exp(-2*core_radius^2/w^2)-1)...
    -(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb)))/comTerm2(pump,Psf,Psb))...
    -Psf*cs_aS*N0*(1-exp(-2*core_radius^2/w^2))-2*h*freqS*dvS*cs_eS*N0*comTerm1(pump,Psf,Psb)/comTerm2(pump,Psf,Psb)...
    *((exp(-2*core_radius^2/w^2)-1)-(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb)))/comTerm2(pump,Psf,Psb))...
    -loss_b*Psf*(1-exp(-2*core_radius^2/w^2));
%backward signal power
dPsbackward_dz = @(z,pump,Psf,Psb) Psb*(cs_eS+cs_aS)*N0*comTerm1(pump,Psf,Psb)/comTerm2(pump,Psf,Psb)*((exp(-2*core_radius^2/w^2)-1)...
    -(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb)))/comTerm2(pump,Psf,Psb))...
    +Psb*cs_aS*N0*(1-exp(-2*core_radius^2/w^2))+2*h*freqS*dvS*cs_eS*N0*comTerm1(pump,Psf,Psb)/comTerm2(pump,Psf,Psb)...
    *((exp(-2*core_radius^2/w^2)-1)-(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb)))/comTerm2(pump,Psf,Psb))...
    +loss_b*Psb*(1-exp(-2*core_radius^2/w^2));

z = 0:dz:len;  
Pp = zeros(1,length(z)); 
Psforward = zeros(1,length(z)); 
Psbackward = zeros(1,length(z)); 

PSout = zeros(1,length(Pp0)); %stores total signal power emitted from fiber for each Pp0
Pabs = zeros(1,length(Pp0)); %stores values for pump power absorbed for each Pp0
dQ_dt1 = zeros(1,length(Pp0)); %stores values for dQ/dt for each Pp0
dQz_dt2 = zeros(length(Pp0),length(z)); %stores dQ/dt as a function of z for each Pp0
dQz_dtSE = zeros(length(Pp0),length(z)); %stores the SE contribution to dQ/dt as a function of z for each Pp0
dQz_dtPhononLoss = zeros(length(Pp0),length(z)); %stores the phonon loss contribution to dQ/dt as a function of z for each Pp0
dQz_dtScat = zeros(length(Pp0),length(z)); %stores the scattering loss contribution to dQ/dt as a function of z for each Pp0
dSF_dz = zeros(length(Pp0),length(z)); %stores the forward signal contribution to dQ/dt as a function of z for each Pp0
dSB_dz = zeros(length(Pp0),length(z)); %stores the backward signal contribution to dQ/dt as a function of z for each Pp0 
dP_dz = zeros(length(Pp0),length(z));%stores the Pump abs contribution to dQ/dt as a function of z for each Pp0
dQ_dt2 = zeros(1,length(Pp0)); %stores values for dQ/dt for each Pp0
dQ_dtSE = zeros(1,length(Pp0)); %stores the SE contribution to dQ/dt for each Pp0
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
dTzSE = zeros(length(Pp0),length(z)); %stores the change in temp at a location z due to SE
dTzP = zeros(length(Pp0),length(z)); %stores the change in temp at a location z due to the pump
dTzSF = zeros(length(Pp0),length(z)); %stores the change in temp at a location z due to the forward sig
dTzSB = zeros(length(Pp0),length(z)); %stores the change in temp at a location z due to the backward sig
dTavg = zeros(1,length(Pp0)); %stores the average change in temperature across the fiber 
dTmin = zeros(1,length(Pp0)); 
dTmax = zeros(1,length(Pp0)); 



for p = 1:length(Pp0)

errorf = NaN(1,maxIterations);
forward = NaN(1,maxIterations);
    
%calculate initial guess for Psforward(1)
nps = Area/gamma/tauRad/(cs_aP+cs_eP);
nss = Area/gamma/tauRad/(cs_aS+cs_eS);
np0 = Pp0(p)/h/freqP;
Psforward(1) = ((gamma*cs_aS*N0*len-log(R1*R2)/2)*nss-np0*(1-exp(-gamma*cs_aP*N0*len)*exp((gamma*cs_aS*N0*len-log(R1*R2)/2)*nss/nps)))...
    /(1-1/sqrt(R1*R2)+(sqrt(R1*R2)-1)/R1)*h*freqS;

csRatio = (cs_aP+cs_eP)/(cs_aS+cs_eS);
minR2 = exp(2*cs_aS*N0*len-2*(cs_aP*N0*len-log(1))/csRatio)
    
Psbackward(1) = Psforward(1)/R1;
Pp(1) = Pp0(p); %W  

iterations = 0;

while iterations <= maxIterations
 
    %forward propagation
    dz = abs(dz);
    for i=1:(length(z)-1) % calculation loop

        RK1p = dPp_dz(z(i),Pp(i),Psforward(i),Psbackward(i));
        RK1f = dPsforward_dz(z(i),Pp(i),Psforward(i),Psbackward(i)); 
        RK1b = dPsbackward_dz(z(i),Pp(i),Psforward(i),Psbackward(i));

        RK2p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(i)+0.5*dz*RK1f,Psbackward(i)+0.5*dz*RK1b);
        RK2f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(i)+0.5*dz*RK1f,Psbackward(i)+0.5*dz*RK1b);
        RK2b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(i)+0.5*dz*RK1f,Psbackward(i)+0.5*dz*RK1b);

        RK3p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(i)+0.5*dz*RK2f,Psbackward(i)+0.5*dz*RK2b);
        RK3f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(i)+0.5*dz*RK2f,Psbackward(i)+0.5*dz*RK2b);
        RK3b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(i)+0.5*dz*RK2f,Psbackward(i)+0.5*dz*RK2b);

        RK4p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(i)+0.5*dz*RK3f,Psbackward(i)+0.5*dz*RK3b);
        RK4f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(i)+0.5*dz*RK3f,Psbackward(i)+0.5*dz*RK3b);
        RK4b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(i)+0.5*dz*RK3f,Psbackward(i)+0.5*dz*RK3b);

        Pp(i+1) = Pp(i) + (1/6)*(RK1p+2*RK2p+2*RK3p+RK4p)*dz;
        Psforward(i+1) = Psforward(i) + (1/6)*(RK1f+2*RK2f+2*RK3f+RK4f)*dz;  % main equation
        Psbackward(i+1) = Psbackward(i) + (1/6)*(RK1b+2*RK2b+2*RK3b+RK4b)*dz;  % main equation
    end
    
    iterations = iterations+1;
    
    %check/modify boundary conditions
    %errorf(iterations) = Psforward(length(z))*R2-Psbackward(length(z));
    errorf(iterations) = Psforward(length(z))*R2/Psbackward(length(z))-1;
    if abs(errorf(iterations)) < maxError || iterations > maxIterations
        break
    else
        forward(iterations) = Psforward(1);
        Psforward(1) = Psforward(1)+Psforward(1)*errorStep*errorf(iterations)/abs(errorf(iterations)); 
        Psbackward(1) = Psforward(1)/R1;
    end
    
    %backward propagation
    dz = -abs(dz);
    for i=flip(2:length(z)) % calculation loop
        
        RK1p = dPp_dz(z(i),Pp(i),Psforward(i),Psbackward(i));
        RK1f = dPsforward_dz(z(i),Pp(i),Psforward(i),Psbackward(i)); 
        RK1b = dPsbackward_dz(z(i),Pp(i),Psforward(i),Psbackward(i));

        RK2p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(i)+0.5*dz*RK1f,Psbackward(i)+0.5*dz*RK1b);
        RK2f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(i)+0.5*dz*RK1f,Psbackward(i)+0.5*dz*RK1b);
        RK2b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(i)+0.5*dz*RK1f,Psbackward(i)+0.5*dz*RK1b);

        RK3p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(i)+0.5*dz*RK2f,Psbackward(i)+0.5*dz*RK2b);
        RK3f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(i)+0.5*dz*RK2f,Psbackward(i)+0.5*dz*RK2b);
        RK3b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(i)+0.5*dz*RK2f,Psbackward(i)+0.5*dz*RK2b);

        RK4p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(i)+0.5*dz*RK3f,Psbackward(i)+0.5*dz*RK3b);
        RK4f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(i)+0.5*dz*RK3f,Psbackward(i)+0.5*dz*RK3b);
        RK4b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(i)+0.5*dz*RK3f,Psbackward(i)+0.5*dz*RK3b);
        
        Pp(i-1) = Pp(i) + (1/6)*(RK1p+2*RK2p+2*RK3p+RK4p)*dz;
        Psforward(i-1) = Psforward(i) + (1/6)*(RK1f+2*RK2f+2*RK3f+RK4f).'*dz;  % main equation
        if sum(Psforward(:,i-1)<zeros(length(lambdaS),1)) > 0
            Psforward(Psforward(:,i-1)<zeros(length(lambdaS),1),i-1) = 0;
        end
        Psbackward(i-1) = Psbackward(i) + (1/6)*(RK1b+2*RK2b+2*RK3b+RK4b).'*dz;  % main equation
        if sum(Psbackward(:,i-1)<zeros(length(lambdaS),1)) > 0
            Psbackward(Psbackward(:,i-1)<zeros(length(lambdaS),1),i-1) = 0;
        end
    end
    
    Pp(1) = Pp0(p);
    Psforward(1) = forward(iterations)+forward(iterations)*errorStep*errorf(iterations)/abs(errorf(iterations)); 
    Psbackward(1) = Psforward(1)/R1;
end

%calculate gain 
gain = zeros(1,length(z));
numPoints = 5;
dr = core_radius/numPoints;
for r = 0:dr:core_radius-dr %takes the average gain across radius of fiber
    for i = 1:length(z)
        gain(i) = gain(i)+gammaS(Pp(i),Psforward(i)+Psbackward(i),r)/numPoints;
    end
end
avgGain = mean(gain); %stores the average gain across the fiber 
maxGain = max(gain); %stores the maximum gain across the fiber

%calculate dQ/dt
%method 2: solving differential equation
T1 = @(Ps) Ps*cs_aS;
T2 = @(Ps) Ps*(cs_aS+cs_eS);
%find dQ/dt as a function of z
for zz = 1:length(z)
    dP_dz(p,zz) = -N0*comTerm1(Pp(zz),Psforward(zz),Psbackward(zz))/comTerm2(Pp(zz),Psforward(zz),Psbackward(zz))...
        *Pp(zz)*(cs_aP+cs_eP)*((exp(-2*core_radius^2/w^2)-1)...
        -(1+tauRatio)*log((tauRatio+1+comTerm2(Pp(zz),Psforward(zz),Psbackward(zz))*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(zz),Psbackward(zz))))...
        /comTerm2(Pp(zz),Psforward(zz),Psbackward(zz)))...
        -N0*Pp(zz)*cs_aP*(1-exp(-2*core_radius^2/w^2));
    
    dSF_dz(p,zz) = -N0*comTerm1(Pp(zz),Psforward(zz),Psbackward(zz))/comTerm2(Pp(zz),Psforward(zz),Psbackward(zz))...
        *(T2(Psforward(zz))+2*h*cs_eS*freqS*dvS)*((exp(-2*core_radius^2/w^2)-1)...
        -(1+tauRatio)*log((tauRatio+1+comTerm2(Pp(zz),Psforward(zz),Psbackward(zz))*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(zz),Psbackward(zz))))...
        /comTerm2(Pp(zz),Psforward(zz),Psbackward(zz)))...
        -N0*T1(Psforward(zz))*(1-exp(-2*core_radius^2/w^2));

    dSB_dz(p,zz) = N0*comTerm1(Pp(zz),Psforward(zz),Psbackward(zz))/comTerm2(Pp(zz),Psforward(zz),Psbackward(zz))...
        *(T2(Psbackward(zz))+2*h*cs_eS*freqS*dvS)*((exp(-2*core_radius^2/w^2)-1)...
        -(1+tauRatio)*log((tauRatio+1+comTerm2(Pp(zz),Psforward(zz),Psbackward(zz))*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(zz),Psbackward(zz))))...
        /comTerm2(Pp(zz),Psforward(zz),Psbackward(zz)))...
        +N0*T1(Psbackward(zz))*(1-exp(-2*core_radius^2/w^2));

    dQz_dtSE(p,zz) = -Am*N0*comTerm1(Pp(zz),Psforward(zz),Psbackward(zz))/comTerm2(Pp(zz),Psforward(zz),Psbackward(zz))*h*freqF_SE/tauRad...
        *log((tauRatio+1+comTerm2(Pp(zz),Psforward(zz),Psbackward(zz))*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(zz),Psbackward(zz))));

    dQz_dtPhononLoss(p,zz) = (Pp(zz)+Psforward(zz)+Psbackward(zz))*loss_ba*(1-exp(-2*core_radius^2/w^2));
    
    dQz_dtScat(p,zz) = (Pp(zz)+Psforward(zz)+Psbackward(zz))*loss_bs*(1-exp(-2*core_radius^2/w^2));
    
    dQz_dt2(p,zz) = -dP_dz(p,zz)-dSF_dz(p,zz)+dSB_dz(p,zz)-dQz_dtSE(p,zz)+dQz_dtPhononLoss(p,zz);
end
[maxdQz_dt2(1,p), maxdQz_dt2(2,p)] = max(dQz_dt2(p,:));
%integrate over z to find dQ/dt for whole fiber
dQ_dt2(p) = trapz(z,dQz_dt2(p,:));
dQ_dtSE(p) = trapz(z,-dQz_dtSE(p,:));
dQ_dtPhononLoss(p) = trapz(z,dQz_dtPhononLoss(p,:));
dQ_dtScat(p) = trapz(z,-dQz_dtScat(p,:));
dQ_dtSF(p) = trapz(z,-dSF_dz(p,:));
dQ_dtSB(p) = trapz(z,dSB_dz(p,:));
dQ_dtP(p) = trapz(z,-dP_dz(p,:));
%method 1: Pabs-Ps 
PSout(p) = Psforward(length(z))*t2
Pabs(p) = Pp(1)-Pp(length(z));
dQ_dt1(p) = Pabs(p)-PSout(p)+dQ_dtSE(p)+dQ_dtScat(p);

%calculate change in temperature
b = 63.5e-6; %m
thermalC = 81.4; %W/m^2/K
dTz(p,:) = dQz_dt2(p,:)/2/pi/b/thermalC; %total temperature change
dTzSE(p,:) = -dQz_dtSE(p,:)/2/pi/b/thermalC; %cooling due to SE
dTzP(p,:) = (-dP_dz(p,:)+dQz_dtPhononLoss(p,:))/2/pi/b/thermalC; %heating due to pump and signal absorption (absorption on impurities and yb)
dTzSF(p,:) = -dSF_dz(p,:)/2/pi/b/thermalC; %cooling due to energy removal through forward signal
dTzSB(p,:) = dSB_dz(p,:)/2/pi/b/thermalC; %cooling due to energy removal through backward signal
dTavg(p) = mean(dTz(p,:)); %K
dTmin(p) = min(dTz(p,:)); %K
dTmax(p) = max(dTz(p,:)); %K

%calculate heat extraction for looped fiber case
%r_looped = determineRadius(loops,b);
%h_looped = 40; %W/m^2/K
loops = 20;
segment = floor(length(z)/loops*100);
dQ_looped = zeros(1,segment);
zNew = linspace(0,len,length(z)*100);
dQz_dtInterp = interp1(z,dQz_dt2(p,:),zNew);
for i = 1:loops
    seg = (i-1)*segment+1:i*segment;
    dQ_looped = dQ_looped+dQz_dtInterp(seg);
end
lastBit = loops*segment+1:length(dQz_dtInterp);
zeroarray = zeros(1,segment-length(lastBit));
dQ_looped = dQ_looped+[dQz_dtInterp(lastBit) zeroarray];
dQcheck = trapz(zNew(1:segment),dQ_looped)
%newTempSeg = dQ_looped/2/pi/r_looped/h_looped;

% dT_looped = [];
% for i = 1:loops+1
%     dT_looped = [dT_looped newTempSeg];
% end
% dT_looped = dT_looped(1:length(dQz_dt2));

%calculate N2(z)
N2z = zeros(1,length(z));
for r = 0:dr:core_radius-dr
    for i = 1:length(z)
        N2z(i) = N2z(i)+N2(Pp(i),Psforward(i)+Psbackward(i),r)/numPoints;
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


%graphs
% Pump power
figure(1)
grid on
hold on
plot(z,Pp);
xlabel('z (m)');
ylabel('Pp(z) (W)');
title('Change in pump Power along fiber');

%pump and signals
figure(2)
hold on 
grid on
plot(z,Pp);
plot(z,Psforward);
plot(z,Psbackward);
xlabel('z (m)');
ylabel('Power (W)');
title('Change in Power along fiber');
legend('forward signal','backward signal');
% 

%Signal Gain
figure(3) 
grid on
hold on
plot(z,gain);
xlabel('z (m)');
ylabel('gain(z) (1/m)');
title('Signal gain along fiber');
% 
%N2(z)
figure(4)
grid on
hold on
plot(z,N2z/N0,'DisplayName',sprintf('%g mW',Pp0(p)*1e3));
%plot(z,N0*ones(length(z)))
xlabel('z (m)');
ylabel('N2(z)/N0');
title('N2 along fiber (85mW pump @ 930nm)');

%dQ/dt as a function of z
figure(5)
grid on 
hold on
plot(z,dQz_dt2(p,:),'DisplayName',sprintf('%g mW',Pp0(p)*1e3));
xlabel('z (m)');
ylabel('dQ/dt (J/s/m)');
title('dQ/dt along the fiber');

%dT as a function of z
figure(6)
hold on
grid on
plot(z,dTz(p,:),'DisplayName',sprintf('%g mW',Pp0(p)*1e3));
xlabel('z (m)');
ylabel('Change in Temperature (K)');
title('Change in Temperature Along the Fiber');

% %dT as a function of z for looped case
% figure(7)
% grid on 
% hold on
% plot(z,dT_looped,'DisplayName',sprintf('%g mW',Pp0(p)*1e3));
% plot(z,dTz(p,:),'DisplayName',sprintf('%g mW',Pp0(p)*1e3));
% xlabel('z (m)');
% ylabel('Change in Temperature (K)');
% title('Change in Temperature along the fiber for looped case');

%dQ as a function of z for loop
figure(7)
grid on 
hold on
plot(zNew(1:segment),dQ_looped,'DisplayName',sprintf('%g mW',Pp0(p)*1e3));
plot(z,dQz_dt2(p,:),'DisplayName',sprintf('%g mW',Pp0(p)*1e3));
xlabel('z (m)');
ylabel('Heat extracted per unit length (W/m)');
title('Heat extraction along the fiber for looped case');

%dQ as a function of z for loop of normalized length
figure(8)
grid on 
hold on
plot(zNew(1:segment)/zNew(segment)*100,dQ_looped,'DisplayName',sprintf('%g mW',Pp0(p)*1e3));
xlabel('normalized length (%)');
ylabel('Heat extracted per unit length (W/m)');
title('Heat extraction along the fiber for looped case');

% %dT due to SE as a function of z
% figure(7)
% hold on
% grid on
% plot(z,dTzSE(p,:),'DisplayName',sprintf('%g mW',Pp0(p)*1e3));
% xlabel('z (m)');
% ylabel('dT (K)');
% title('change in Temperature due to SE along the fiber');

% %contributions to dT
% figure(8)
% hold on
% plot(z,dTzSE(p,:))
% plot(z,dTzP(p,:))
% plot(z,dTzP(p,:)+dTzSF(p,:)+dTzSB(p,:))
% plot(z,dTzSF(p,:)+dTzSB(p,:))
% plot(z,dTz(p,:),'lineWidth',3)
% xlabel('z (m)');
% ylabel('Change In Temperature (K)');
% legend('SE','Pump','Pump-Signal','Signal','total');
% title('Contributions to Change in Temperature');
% grid on

% %contributions to dT
% figure(10)
% hold on
% plot(z,dTzSE(p,:))
% plot(z,dTzP(p,:)+dTzSF(p,:)+dTzSB(p,:))
% plot(z,dTz(p,:),'lineWidth',3)
% xlabel('z (m)');
% ylabel('Change In Temperature (K)');
% legend('SE','Absorbed Pump - Signal','total');
% title('Contributions to Change in Temperature');
% grid on

maxSE = min(dTzSE(p,:))

%convergence error
figure(12)
hold on
plot(1:length(isfinite(errorf)),errorf(1:length(isfinite(errorf))));
xlabel('iteration');
ylabel('error')
title('error vs. iteration')
grid on

%change in forward signal with convergence
figure(13)
hold on
plot(1:length(isfinite(forward)),forward(1:length(isfinite(forward))));
xlabel('iteration');
ylabel('forward signal')
title('forward signal vs. iteration')
% grid on

end


% figure(4) %N2(z)
% legend('show');
% figure(6) %dTz
% legend('show');
% figure(5) %dQ/dt(z)
% legend('show')
% plot(z(maxdQz_dt2(2,:)),maxdQz_dt2(1,:));

% %total signal power AND max temperature vs. absorped pump power
% figure(9)
% hold on
% yyaxis left
% plot(Pabs,PSout)
% ylabel('Output Signal (W)');
% yyaxis right
% plot(Pabs,dTmax);
% ylabel('Maximum temperature (K)');
% xlabel('Absorbed Pump Power (W)');
% title('Output Signal vs. Absorbed Pump Power');
% grid on
% 
% %total signal power AND avg temperature vs. Pp0
% figure(10)
% hold on
% yyaxis left
% plot(Pp0,PSout)
% ylabel('Output Signal (W)');
% yyaxis right
% plot(Pp0,dTavg);
% ylabel('Average temperature (K)');
% xlabel('Input Pump Power (W)');
% title('Output Signal vs. Input Pump Power');
% grid on
% 
% %total signal power AND total heat generation vs. Pp0
% figure(11)
% hold on
% yyaxis left
% plot(Pp0,PSout)
% ylabel('Output Signal (W)');
% yyaxis right
% plot(Pp0,dQ_dt2*1e3);
% ylabel('total heat generation (mW)');
% xlabel('Input Pump Power (W)');
% title('Output Signal vs. Input Pump Power');
% grid on
% 
% %total heat generation vs. Pout
% figure(12)
% hold on
% plot(PSout,dQ_dt2*1e3)
% xlabel('Laser output power (W)');
% ylabel('Total heat generation (mW)');
% box on

%maximum heating vs. absorbed pump power
%%total signal power vs. absorped pump power
% figure(10)
% yaxis = plotyy(Pabs*1e3,dTmax,Pabs*1e3,PSout*1e3);
% dTmaxTicks = round(linspace(dTmax(1),dTmax(end),10),2);
% set(yaxis(1),'YTick',dTmaxTicks,'ylim',[dTmaxTicks(1) dTmaxTicks(end)],'Box','off');
% PSoutYaxis = ylim(yaxis(2));
% set(yaxis(2),'YTick',round(linspace(PSoutYaxis(1),PSoutYaxis(end),10),1));
% ylabel(yaxis(1),'change in Temp (K)');
% ylabel(yaxis(2),'Output Signal (mW)');
% xlabel('Pabs (mW)');
% title('Max heating and laser power vs. absorbed pump power (pump = 930nm)');
% grid on

% 
% %Q vs. Pp0
% figure(11)
% grid on
% hold on
% %plot(Pp0,dQ_dt1);
% plot(Pp0,dQ_dt2);
% xlabel('Pp(0) (W)');
% ylabel('dQ/dt (J/s)')
% title('Q vs Initial pump power');
% %legend('Pabs-Ps','differential equation');

% %SE and ASE contribution to Q
% figure(12)
% grid on
% hold on
% plot(Pp0,Ptotsignal);
% %plot(Pp0,dQ_dtSE);
% xlabel('Pp(0) (W)');
% ylabel('total Signal (W)')
% title('Total signal vs. Initial pump power');
% %legend('ASE','other SE');



expectedSlope = t2/(2*loss_b*len+t1+t2)*lambdaP/lambdaS

% %total signal power vs. absorped pump power
% figure(10)
% hold on
% absPow = [10 80];
% plot(absPow,r299.coeff(1)*absPow+r299.coeff(2))
% plot(absPow,r250.coeff(1)*absPow+r250.coeff(2))
% xlabel('Absorbed Pump Power (mW)');
% ylabel('Output Signal (mW)');
% title('Output Signal vs. Absorbed Pump Power');
% grid on
% legend('R2 = 99%','R2 = 50%')



% me = (Pp(length(z))-Pp(1))/len(1);
% ma = (approx(length(z))-approx(1))/len(1);
% mratio = me/ma
% logFactor = (1+exp(-2*core_radius^2/w^2))/2
% slopefactor = log(logFactor)


% sum1 = (exp(-c*h*E21/kT)+exp(-c*h*E22/kT)+exp(-c*h*E23/kT));
% sum2 = exp(-h*c*E21/kT)*(E21-E11)+exp(-h*c*E21/kT)*(E21-E12)+exp(-h*c*E21/kT)*(E21-E13)+exp(-h*c*E21/kT)*(E21-E14)...
%     +exp(-h*c*E22/kT)*(E22-E11)+exp(-h*c*E22/kT)*(E22-E12)+exp(-h*c*E22/kT)*(E22-E13)+exp(-h*c*E22/kT)*(E22-E14)...
%     +exp(-h*c*E23/kT)*(E23-E11)+exp(-h*c*E23/kT)*(E23-E12)+exp(-h*c*E23/kT)*(E23-E13)+exp(-h*c*E23/kT)*(E23-E14);
% lambdatest1 = g1*sum1/sum2;
% 
% lambdatest2 = N2z(11)/(N2j(1,11)*(E21-E11)+N2j(1,11)*(E21-E12)+N2j(1,11)*(E21-E13)+N2j(1,11)*(E21-E14)...
%     +N2j(2,11)*(E22-E11)+N2j(2,11)*(E22-E12)+N2j(2,11)*(E22-E13)+N2j(2,11)*(E22-E14)...
%     +N2j(3,11)*(E23-E11)+N2j(3,11)*(E23-E12)+N2j(3,11)*(E23-E13)+N2j(3,11)*(E23-E14));
% 
% for zz = 1:length(z)
%     sum1 = (E21*exp(-c*h*E21/kT)+E22*exp(-c*h*E22/kT)+E23*exp(-c*h*E23/kT)+(E21-E12)*exp(-c*h*E21/kT)...
%         +(E22-E12)*exp(-c*h*E22/kT)+(E23-E12)*exp(-c*h*E23/kT)+(E21-E13)*exp(-c*h*E21/kT)+(E22-E13)*exp(-c*h*E22/kT)...
%         +(E23-E13)*exp(-c*h*E23/kT)+(E21-E14)*exp(-c*h*E21/kT)+(E22-E14)*exp(-c*h*E22/kT)+(E23-E14)*exp(-c*h*E23/kT));
%     sum2 = exp(-c*h*E21/kT)+exp(-c*h*E22/kT)+exp(-c*h*E23/kT);
%     lambdaFz_dig(zz) = sum2/sum1;
% end
% lambdaF_dig(p) = mean(lambdaFz_dig);


