%numerical result for change in power across fiber for one pump and
%broadboand signal. Automatically converges. Does NOT approximate intensity in
%fiber as P(z)/Area.
%Automatically discards signal wavelengths for which gain is too low
%includes non-radiative term and loss

%close all;
clear all;

%set(0,'DefaultAxesFontSize',18)

len = 20;%.01:0.03:.5;
numdz = 400;
dz = len/numdz;
maxIterations = 5;
maxTries = 1;
gainCutOff = 0.08; %if a signal wavelength has a gain factor that is less 
                   %than gainCutOff percent of the max gain factor, then
                   %the signal at that wavelength is set to 0;
minGain = -2500; %if a signal wavelength has a gain that is less 
               %than minGain than the signal at that wavelength is set to 0
maxError = 1e-8;
eliminationFactor = 1e4;

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

core_radius = 3.1e-6; 
Area = pi*core_radius^2; %m^2, pump area
tauRad = 1.7e-3;
c = 3e8;
h = 6.63e-34; %J*s
kT = 4.11e-21; %J

%cross sectional areas 
wavelengths = (870:1100)*1e-9;
cs_absRAW = xlsread('abs_ZBLANP_MC.xlsx');
cs_abs = [wavelengths; interp1(cs_absRAW(:,1),cs_absRAW(:,2),wavelengths)].';
cs_abs(isnan(cs_abs)) = 0;
cs_emsRAW = xlsread('emm_ZBLANP_MC.xlsx');
cs_ems = [wavelengths; interp1(cs_emsRAW(:,1),cs_emsRAW(:,2),wavelengths)].';
cs_ems(isnan(cs_ems)) = 0;

%calculate mean flourecent wavelength using Mina's formula
integral1 = cs_emsRAW(:,2)./cs_emsRAW(:,1).^4;%cs_emsRAW(:,1).*cs_emsRAW(:,2)./cs_emsRAW(:,1).^5;
integral2 = cs_emsRAW(:,2)./cs_emsRAW(:,1).^5;
lambdaF_SE = trapz(cs_emsRAW(:,1),integral1)/trapz(cs_emsRAW(:,1),integral2);
freqF_SE = c/lambdaF_SE;

%create a variable to help me compare cross section values
crossSection = cs_abs;
crossSection(:,3) = cs_ems(:,2);
crossSection(:,4) = cs_abs(:,2)./cs_ems(:,2);

%waveleghth of laser/cooling pump
dlam = 5e-9; %resolution of the signal spectrum
lambdaP = 1030e-9;%1025e-9:5e-9:1050e-9;%980e-9:2e-9:1050e-9;%860e-9:10e-9:1050e-9;
lambdaS = [1030e-9:dlam:1100e-9]; 
lambdaGains = [1030e-9 1035e-9 1050e-9];%[round(lambdaF_SE*1e9)/1e9 1020e-9 1030e-9 1050e-9];
indexGain = NaN(size(lambdaGains));
for gainLam = 1:length(lambdaGains)
    indexGain(gainLam) = find(round(lambdaS*1e9 - lambdaGains(gainLam)*1e9) == 0);
end
%validWavelengths = find(ones(1,length(lambdaS)));
freqP = c./lambdaP;
freqS = c./lambdaS;
%dvS = zeros(1,length(lambdaS));
dvS = zeros(1,length(lambdaS));%c./(lambdaS-dlam/2)-c./(lambdaS+dlam/2); %zeros(1,length(lambdaS));%
%cross sectional areas for a given wavelength

indexP = NaN(size(lambdaP));
for pumpLam = 1:length(lambdaP)
    indexP(pumpLam) = find(round(cs_abs(:,1)*1e9 - lambdaP(pumpLam)*1e9) == 0);
end
cs_aP = cs_abs(indexP,2).';
cs_eP = cs_ems(indexP,2).';

indexS = NaN(size(lambdaS));
for i = 1:length(lambdaS)
    indexS(i) = find(round(cs_abs(:,1)*1e9 - lambdaS(i)*1e9) == 0);
end
cs_aS = cs_abs(indexS,2).';
cs_eS = cs_ems(indexS,2).';


lam = 1000e-9;
NA = .13;
%V = 2*pi*core_radius/lam*NA;
V = 2.44;
w = core_radius*(0.65+1.619/V^1.5+2.879/V^6);
etta = 1-exp(-2*core_radius^2/w^2);


N0_wtPercent = 1;%[0.00005:.00005:.001 .001:.005:.1 .1:.1:3 3:15];
N0 = N0_wtPercent*2.42e26;  %m^-3
Nc = 6.47e27; %m^-3
tauNonRad = 2*pi/9*tauRad*(Nc/N0)^2;%67.6e-3;%211.65e-3;%1e10;%
tauRatio = tauRad/tauNonRad;
tau = tauRad*tauNonRad/(tauRad+tauNonRad);
f = 1/Area;
gamma = 1;
totLoss = .02; %dB/m
loss_b = totLoss/etta*log(10)/10; %m^-1
loss_bs = 1*loss_b/4;%m^-1
loss_ba = 3*loss_b/4;%0; %m^-1
Am = pi*w^2/2;


z = 0:dz:len;  
Pp = zeros(1,length(z)); 
Psforward = zeros(length(lambdaS),length(z)); 
Psbackward = zeros(length(lambdaS),length(z)); 

Pp0 = zeros(1,length(N0)); %stores optimal Pp0 for each wavelength
Popt = zeros(1,length(N0)); %stores optimal Pp0 for each wavelength
Ptotsignal = zeros(1,length(lambdaP)); %stores total signal power emitted from fiber for each lambdaP
Pabs = zeros(1,length(lambdaP)); %stores values for pump power absorbed for each lambdaP
dQ_dt1 = zeros(1,length(lambdaP)); %stores values for dQ/dt for each lambdaP
dQz_dt2 = zeros(length(lambdaP),length(z)); %stores dQ/dt as a function of z for each lambdaP
dQz_dtSE = zeros(length(lambdaP),length(z)); %stores the SE contribution to dQ/dt as a function of z for each lambdaP
dQzcorr_dtSE = zeros(length(lambdaP),length(z)); %stores the correction term the SE cooling as a funciton of z for each Pp0
dQz_dtPhononLoss = zeros(length(lambdaP),length(z)); %stores the phonon loss contribution to dQ/dt as a function of z for each lambdaP
dQz_dtScat = zeros(length(lambdaP),length(z)); %stores the scattering loss contribution to dQ/dt as a function of z for each lambdaP
dSF_dz = zeros(length(lambdaP),length(z)); %stores the forward signal contribution to dQ/dt as a function of z for each lambdaP
dSB_dz = zeros(length(lambdaP),length(z)); %stores the backward signal contribution to dQ/dt as a function of z for each lambdaP 
dP_dz = zeros(length(lambdaP),length(z));%stores the Pump abs contribution to dQ/dt as a function of z for each lambdaP
dQ_dt2 = zeros(1,length(lambdaP)); %stores values for dQ/dt for each lambdaP
dQ_dtSE = zeros(1,length(lambdaP)); %stores the SE contribution to dQ/dt for each lambdaP
dQ_dtSEcorr = zeros(1,length(lambdaP)); %stores the SE correction term to dQ/dt for each Pp0
dQ_dtPhononLoss = zeros(1,length(lambdaP)); %stores the Phonon loss contribution to dQ/dt for each lambdaP
dQ_dtScat = zeros(1,length(lambdaP)); %stores the scattering loss contribution to dQ/dt for each lambdaP
dQ_dtSF = zeros(1,length(lambdaP)); %stores the forward signal contribution to dQ/dt for each lambdaP
dQ_dtSB = zeros(1,length(lambdaP)); %stores the backward signal contribution to dQ/dt for each lambdaP
dQ_dtP = zeros(1,length(lambdaP)); %stores the pump abs contribution to dQ/dt for each lambdaP
avgGain = zeros(length(lambdaP),length(lambdaS)); %stores the average gain across the fiber for each wavelength
maxGain = zeros(length(lambdaP),length(lambdaS)); %stores the maximum gain across the fiber for each wavelength
maxdQz_dt2 = zeros(2,length(lambdaP)); %stores the coordinates of the maximum value of dQz_dt2
                                    %first row: max values
                                    %second row: index of z where max value is
avgdQz_dt2 = zeros(1,length(lambdaP));                                    
lambdaF_SL = zeros(1,length(lambdaP)); %stores the mean flouresnce wavelength for each pump power calculated from the 
                                    %populations of each sublevel.
lambdaFF_ASE = zeros(1,length(lambdaP)); %stores the mean flouresence wavelength of the ASE for each pump wavelength
lambdaFB_ASE = zeros(1,length(lambdaP));
dTz = zeros(length(lambdaP),length(z));  %stores the change in temp at a location z along the fiber
dTavg = zeros(1,length(lambdaP)); %stores the average change in temperature across the fiber 
dTmin = zeros(1,length(lambdaP)); %stores the minimum change in temperature across the fiber (most cooling if negative)

for l = 1:length(lambdaP)
    
IsatPa = h*freqP(l)/cs_aP(l)/tauRad;
PsatPa = IsatPa*pi*w^2/2;
IsatSa = h*freqS./cs_aS/tauRad;
IsatP = h*freqP(l)/(cs_aP(l)+cs_eP(l))/tauRad;
PsatP = IsatP*pi*w^2/2;
IsatS = h*freqS./(cs_aS+cs_eS)/tauRad;

validWavelengths = find(ones(1,length(lambdaS)));

%Gaussian mode intensity function
fr = @(r) 2/pi/w^2*exp(-2*r.^2/w^2);
%N2(z) - Pstot is 1xlength(lambdaS) array where each element is the sum of
%the backward and forward power for that wavelength at position 
N2 = @(Pp,Pstot,r,vw) (Pp*fr(r)/IsatPa+sum(Pstot*fr(r)./IsatSa(vw)))/(tauRatio+1+Pp*fr(r)/IsatP+sum(Pstot*fr(r)./IsatS(vw)))*N0; 
%signal gain coefficients
gammaS = @(Pp,Pstot,r,vw) N2(Pp,Pstot,r,vw)*(cs_eS(vw)+cs_aS(vw))-cs_aS(vw)*N0; %returns and array of size 1xlength(validWavelengths) 
gammaSE = @(Pp,Pstot,r,vw) cs_eS(vw)*N2(Pp,Pstot,r,vw); %returns and array of size 1xlength(validWavelengths)

comTerm1 = @(pump,Psf,Psb,vw) (pump/PsatPa+sum((Psf+Psb)./IsatSa(vw))/Am);
comTerm2 = @(pump,Psf,Psb,vw) (pump/PsatP+sum((Psf+Psb)./IsatS(vw))/Am);
%pump power - Psf and Psb are 1xlength(validWavelengths) arrays where each element
%is the forward and backward power (respectively) for that wavelegth and
%postion.
dPp_dz = @(z,pump,Psf,Psb,vw) pump*cs_aP(l)*N0/comTerm2(pump,Psf,Psb,vw)...
    *(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb,vw)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb,vw)))...
    +pump*N0*(cs_aP(l)*sum((Psf+Psb)./IsatS(vw))/Am-(cs_aP(l)+cs_eP(l))*sum((Psf+Psb)./IsatSa(vw))/Am)/comTerm2(pump,Psf,Psb,vw)...
    *((exp(-2*core_radius^2/w^2)-1)-(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb,vw)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb,vw)))/comTerm2(pump,Psf,Psb,vw))...
    -loss_b*pump*(1-exp(-2*core_radius^2/w^2));
%forward signal power
dPsforward_dz = @(z,pump,Psf,Psb,vw,dv) -Psf.*(cs_eS(vw)+cs_aS(vw))*N0*comTerm1(pump,Psf,Psb,vw)/comTerm2(pump,Psf,Psb,vw)*((exp(-2*core_radius^2/w^2)-1)...
    -(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb,vw)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb,vw)))/comTerm2(pump,Psf,Psb,vw))...
    -Psf.*cs_aS(vw)*N0*(1-exp(-2*core_radius^2/w^2))-2*h*freqS(vw).*dv(vw).*cs_eS(vw)*N0*comTerm1(pump,Psf,Psb,vw)/comTerm2(pump,Psf,Psb,vw)...
    *((exp(-2*core_radius^2/w^2)-1)-(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb,vw)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb,vw)))/comTerm2(pump,Psf,Psb,vw))...
    -loss_b*Psf*(1-exp(-2*core_radius^2/w^2));
%backward signal power
dPsbackward_dz = @(z,pump,Psf,Psb,vw,dv) Psb.*(cs_eS(vw)+cs_aS(vw))*N0*comTerm1(pump,Psf,Psb,vw)/comTerm2(pump,Psf,Psb,vw)*((exp(-2*core_radius^2/w^2)-1)...
    -(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb,vw)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb,vw)))/comTerm2(pump,Psf,Psb,vw))...
    +Psb.*cs_aS(vw)*N0*(1-exp(-2*core_radius^2/w^2))+2*h*freqS(vw).*dv(vw).*cs_eS(vw)*N0*comTerm1(pump,Psf,Psb,vw)/comTerm2(pump,Psf,Psb,vw)...
    *((exp(-2*core_radius^2/w^2)-1)-(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb,vw)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb,vw)))/comTerm2(pump,Psf,Psb,vw))...
    +loss_b*Psb*(1-exp(-2*core_radius^2/w^2));
  

%find pump power such that the fiber has zero net temp change
noTempChange = @(pump) N0*PsatP/PsatPa*log((tauRatio+1+pump/PsatP*exp(-2*core_radius^2/w^2))/(tauRatio+1+pump/PsatP))*(Am*h*freqF_SE/tauRad-(cs_aP(l)+cs_eP(l))*PsatP*(1+tauRatio))...
    +pump*etta*(loss_ba+N0*cs_aP(l)-N0*PsatP/PsatPa*(cs_aP(l)+cs_eP(l)));

Popt(l) = IsatP*pi/2*(-2*core_radius^2/log(1-etta))*(-1*(2-etta)+sqrt(etta^2+4*(1-etta)*cs_aP(l)*N0/loss_ba*(tau/tauRad*lambdaP(l)/lambdaF_SE-1)))/2/(1-etta);

Pp0(l) = fzero(noTempChange,Popt(l)*10);
if isnan(Pp0(l))
    Pp0(l) = fzero(noTempChange,Popt(l)*20);
end

Pp(1) = Pp0(l); %W  

iterations = 0;
tries = 0;


while iterations <= maxIterations
 
    %forward propagation
    dz = abs(dz);
    for i=1:(length(z)-1) % calculation loop

        RK1p = dPp_dz(z(i),Pp(i),Psforward(validWavelengths,i).',Psbackward(validWavelengths,i).',validWavelengths);
        RK1f = dPsforward_dz(z(i),Pp(i),Psforward(validWavelengths,i).',Psbackward(validWavelengths,i).',validWavelengths,dvS); %array of 1xlength(lambdaS)
        RK1b = dPsbackward_dz(z(i),Pp(i),Psforward(validWavelengths,i).',Psbackward(validWavelengths,i).',validWavelengths,dvS); %array of 1xlength(lambdaS)

        RK2p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(validWavelengths,i).'+0.5*dz*RK1f,Psbackward(validWavelengths,i).'+0.5*dz*RK1b,validWavelengths);
        RK2f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(validWavelengths,i).'+0.5*dz*RK1f,Psbackward(validWavelengths,i).'+0.5*dz*RK1b,validWavelengths,dvS);
        RK2b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(validWavelengths,i).'+0.5*dz*RK1f,Psbackward(validWavelengths,i).'+0.5*dz*RK1b,validWavelengths,dvS);

        RK3p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(validWavelengths,i).'+0.5*dz*RK2f,Psbackward(validWavelengths,i).'+0.5*dz*RK2b,validWavelengths);
        RK3f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(validWavelengths,i).'+0.5*dz*RK2f,Psbackward(validWavelengths,i).'+0.5*dz*RK2b,validWavelengths,dvS);
        RK3b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(validWavelengths,i).'+0.5*dz*RK2f,Psbackward(validWavelengths,i).'+0.5*dz*RK2b,validWavelengths,dvS);

        RK4p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(validWavelengths,i).'+0.5*dz*RK3f,Psbackward(validWavelengths,i).'+0.5*dz*RK3b,validWavelengths);
        RK4f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(validWavelengths,i).'+0.5*dz*RK3f,Psbackward(validWavelengths,i).'+0.5*dz*RK3b,validWavelengths,dvS);
        RK4b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(validWavelengths,i).'+0.5*dz*RK3f,Psbackward(validWavelengths,i).'+0.5*dz*RK3b,validWavelengths,dvS);

        Pp(i+1) = Pp(i) + (1/6)*(RK1p+2*RK2p+2*RK3p+RK4p)*dz;
        Psforward(validWavelengths,i+1) = Psforward(validWavelengths,i) + (1/6)*(RK1f+2*RK2f+2*RK3f+RK4f).'*dz;  % main equation
        if sum(Psforward(:,i+1)<zeros(length(lambdaS),1)) > 0
            Psforward(Psforward(:,i+1)<zeros(length(lambdaS),1),i+1) = 0;
        end
        Psbackward(validWavelengths,i+1) = Psbackward(validWavelengths,i) + (1/6)*(RK1b+2*RK2b+2*RK3b+RK4b).'*dz;  % main equation
        if sum(Psbackward(:,i+1)<zeros(length(lambdaS),1)) > 0
            Psbackward(Psbackward(:,i+1)<zeros(length(lambdaS),1),i+1) = 0;
        end
    end
    
    %check/modify boundary conditions
    if sum(abs(Psbackward(validWavelengths,length(z))) > maxError) == 0 && (sum(Psbackward(validWavelengths,length(z))) > 0 || sum(dvS)==0)
        break
    else
        Psbackward(validWavelengths,length(z)) = 0;
    end
    
    %backward propagation
    dz = -abs(dz);
    for i=flip(2:length(z)) % calculation loop
        
        RK1p = dPp_dz(z(i),Pp(i),Psforward(validWavelengths,i).',Psbackward(validWavelengths,i).',validWavelengths);
        RK1f = dPsforward_dz(z(i),Pp(i),Psforward(validWavelengths,i).',Psbackward(validWavelengths,i).',validWavelengths,dvS); %array of 1xlength(lambdaS)
        RK1b = dPsbackward_dz(z(i),Pp(i),Psforward(validWavelengths,i).',Psbackward(validWavelengths,i).',validWavelengths,dvS); %array of 1xlength(lambdaS)

        RK2p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(validWavelengths,i).'+0.5*dz*RK1f,Psbackward(validWavelengths,i).'+0.5*dz*RK1b,validWavelengths);
        RK2f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(validWavelengths,i).'+0.5*dz*RK1f,Psbackward(validWavelengths,i).'+0.5*dz*RK1b,validWavelengths,dvS);
        RK2b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(validWavelengths,i).'+0.5*dz*RK1f,Psbackward(validWavelengths,i).'+0.5*dz*RK1b,validWavelengths,dvS);

        RK3p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(validWavelengths,i).'+0.5*dz*RK2f,Psbackward(validWavelengths,i).'+0.5*dz*RK2b,validWavelengths);
        RK3f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(validWavelengths,i).'+0.5*dz*RK2f,Psbackward(validWavelengths,i).'+0.5*dz*RK2b,validWavelengths,dvS);
        RK3b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(validWavelengths,i).'+0.5*dz*RK2f,Psbackward(validWavelengths,i).'+0.5*dz*RK2b,validWavelengths,dvS);

        RK4p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(validWavelengths,i).'+0.5*dz*RK3f,Psbackward(validWavelengths,i).'+0.5*dz*RK3b,validWavelengths);
        RK4f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(validWavelengths,i).'+0.5*dz*RK3f,Psbackward(validWavelengths,i).'+0.5*dz*RK3b,validWavelengths,dvS);
        RK4b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(validWavelengths,i).'+0.5*dz*RK3f,Psbackward(validWavelengths,i).'+0.5*dz*RK3b,validWavelengths,dvS);

        Pp(i-1) = Pp(i) + (1/6)*(RK1p+2*RK2p+2*RK3p+RK4p)*dz;
        Psforward(validWavelengths,i-1) = Psforward(validWavelengths,i) + (1/6)*(RK1f+2*RK2f+2*RK3f+RK4f).'*dz;  % main equation
        if sum(Psforward(:,i-1)<zeros(length(lambdaS),1)) > 0
            Psforward(Psforward(:,i-1)<zeros(length(lambdaS),1),i-1) = 0;
        end
        Psbackward(validWavelengths,i-1) = Psbackward(validWavelengths,i) + (1/6)*(RK1b+2*RK2b+2*RK3b+RK4b).'*dz;  % main equation
        if sum(Psbackward(:,i-1)<zeros(length(lambdaS),1)) > 0
            Psbackward(Psbackward(:,i-1)<zeros(length(lambdaS),1),i-1) = 0;
        end 
    end
    
    iterations = iterations+1;
    
    %check/modify boundary conditions
    if (sum(abs(Psforward(validWavelengths,1)) > maxError) == 0 && abs(Pp(1)-Pp0) < maxError) %%|| iterations > maxIterations
        break
    end
    
    %if can't converge with all wavelengths, remove troublesome wavelengths
    %that have sufficiently small gain.
    if iterations > maxIterations
        %change Psforward/Psbackward so that data for invalid wavelengths are set to zero.
        %data for valid wavelengths is kept the same.
        invalidWavelengths = ~ismember(1:length(lambdaS),validWavelengths);
        Psforward(invalidWavelengths,:) = zeros(sum(invalidWavelengths),length(z));
        Psbackward(invalidWavelengths,:) = zeros(sum(invalidWavelengths),length(z));
        
        tries = tries+1;
        if tries > maxTries
            break
        end
        
        newValidWavelengths = find(abs(Psforward(validWavelengths,1)) < eliminationFactor*maxError).';
        if isequal(newValidWavelengths,validWavelengths)
            break
        else
            validWavelengths = newValidWavelengths;
            %calculate new dvS
            newdlam = lambdaS(validWavelengths(2)) - lambdaS(validWavelengths(1));
            dvS(1) = c./(lambdaS(validWavelengths(1))-newdlam/2)-c./(lambdaS(validWavelengths(1))+newdlam/2); 
            for lam = 2:length(validWavelengths)-1
                indexB = validWavelengths(lam-1);
                index = validWavelengths(lam);
                indexF = validWavelengths(lam+1);
                newdlamF = lambdaS(indexF)-lambdaS(index);
                newdlamB = lambdaS(index)-lambdaS(indexB);
                dvS(index) = c./(lambdaS(index)-newdlamB/2)-c./(lambdaS(index)+newdlamF/2); 
            end
            newdlam = lambdaS(validWavelengths(end)) - lambdaS(validWavelengths(end-1));
            dvS(validWavelengths(end)) = c./(lambdaS(validWavelengths(end))-newdlam/2)-c./(lambdaS(end)+newdlam/2); 
        end
        iterations = 0;
    end

    Psforward(validWavelengths,1) = 0;
    Pp(1) = Pp0;
end

%change Psforward/Psbackward so that data for invalid wavelengths are set to zero.
%data for valid wavelengths is kept the same.
invalidWavelengths = ~ismember(1:length(lambdaS),validWavelengths);
Psforward(invalidWavelengths,:) = zeros(sum(invalidWavelengths),length(z));
Psbackward(invalidWavelengths,:) = zeros(sum(invalidWavelengths),length(z));

%calculate gain 
gain = zeros(length(lambdaS),length(z));
numPoints = 5;
dr = core_radius/numPoints;
for r = 0:dr:core_radius-dr %takes the average gain across radius of fiber
    for i = 1:length(z)
        gain(:,i) = gain(:,i)+gammaS(Pp(i),Psforward(:,i).'+Psbackward(:,i).',r,1:length(lambdaS)).'/numPoints;
    end
end
for i = 1:length(lambdaS)
    avgGain(l,i) = mean(gain(i,:));
    maxGain(l,i) = max(gain(i,:));
end

%calculate mean flourescent wavelength of ASE spectrum
numWaves = length(lambdaS);
lambdaFF_ASE(l) = 0;
lambdaFB_ASE(l) = 0;
for i = 1:numWaves
    lambdaFF_ASE(l) = lambdaFF_ASE(l)+lambdaS(i)*Psforward(i,length(z));
    lambdaFB_ASE(l) = lambdaFB_ASE(l)+lambdaS(i)*Psbackward(i,1);
end
lambdaFF_ASE(l) = lambdaFF_ASE(l)/sum(Psforward(:,length(z)));
lambdaFB_ASE(l) = lambdaFB_ASE(l)/sum(Psbackward(:,1));

%calculate dQ/dt
%method 2: solving differential equation
T1 = @(Ps) sum(Ps.*cs_aS);
T2 = @(Ps) sum(Ps.*(cs_aS+cs_eS));
E = sum(cs_eS.*freqS.*dvS);
%find dQ/dt as a function of z
awl = 1:length(lambdaS); %indexes of all wavelengths
for zz = 1:length(z)
    dP_dz(l,zz) = -N0*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)...
        *Pp(zz)*(cs_aP(l)+cs_eP(l))*((exp(-2*core_radius^2/w^2)-1)...
        -(1+tauRatio)*log((tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)))...
        /comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl))...
        -N0*Pp(zz)*cs_aP(l)*(1-exp(-2*core_radius^2/w^2));
    
    dSF_dz(l,zz) = -N0*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)...
        *(T2(Psforward(:,zz).')+2*h*sum(cs_eS(validWavelengths).*freqS(validWavelengths).*dvS(validWavelengths)))*((exp(-2*core_radius^2/w^2)-1)...
        -(1+tauRatio)*log((tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)))...
        /comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl))...
        -N0*T1(Psforward(:,zz).')*(1-exp(-2*core_radius^2/w^2));

    dSB_dz(l,zz) = N0*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)...
        *(T2(Psbackward(:,zz).')+2*h*sum(cs_eS(validWavelengths).*freqS(validWavelengths).*dvS(validWavelengths)))*((exp(-2*core_radius^2/w^2)-1)...
        -(1+tauRatio)*log((tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)))...
        /comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl))...
        +N0*T1(Psbackward(:,zz).')*(1-exp(-2*core_radius^2/w^2));

    dQz_dtSE(l,zz) = -Am*N0*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)*h*freqF_SE/tauRad...
        *log((tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)));

    
    dQzcorr_dtSE(l,zz) = -4*h*N0*E*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)...
        *((exp(-2*core_radius^2/w^2)-1)...
        -(1+tauRatio)*log((tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)))...
        /comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl));
    
    dQz_dtPhononLoss(l,zz) = (Pp(zz)+sum(Psforward(:,zz)+Psbackward(:,zz)))*loss_ba*(1-exp(-2*core_radius^2/w^2));
    
    dQz_dtScat(l,zz) = (Pp(zz)+sum(Psforward(:,zz)+Psbackward(:,zz)))*loss_bs*(1-exp(-2*core_radius^2/w^2));
    
    dQz_dt2(l,zz) = -dP_dz(l,zz)-dSF_dz(l,zz)+dSB_dz(l,zz)-dQz_dtSE(l,zz)+dQzcorr_dtSE(l,zz)+dQz_dtPhononLoss(l,zz);
end
[maxdQz_dt2(1,l), maxdQz_dt2(2,l)] = min(dQz_dt2(l,:));
avgdQz_dt2(l) = mean(dQz_dt2(l,:));
%integrate over z to find dQ/dt for whole fiber
dQ_dt2(l) = trapz(z,dQz_dt2(l,:));
dQ_dtSE(l) = trapz(z,-dQz_dtSE(l,:));
dQ_dtSEcorr(l) = trapz(z,dQzcorr_dtSE(l,:));
dQ_dtPhononLoss(l) = trapz(z,dQz_dtPhononLoss(l,:));
dQ_dtScat(l) = trapz(z,-dQz_dtScat(l,:));
dQ_dtSF(l) = trapz(z,-dSF_dz(l,:));
dQ_dtSB(l) = trapz(z,dSB_dz(l,:));
dQ_dtP(l) = trapz(z,-dP_dz(l,:));
%method 1: Pabs-Ps 
Ptotsignal(l) = sum(Psforward(:,length(z))-Psforward(:,1)+Psbackward(:,1)-Psbackward(:,length(z)));
Pabs(l) = Pp(1)-Pp(length(z));
dQ_dt1(l) = Pabs(l)-Ptotsignal(l)+dQ_dtSE(l)+dQ_dtSEcorr(l)+dQ_dtScat(l);

%calculate change in temperature
b = 67.3e-6; %m
thermalC = 81.4; %W/m^2/K (for a cladding radius of 62.5 um)
dTz(l,:) = dQz_dt2(l,:)/2/pi/b/thermalC; %K
dTavg(l) = mean(dTz(l,:)); %K
dTmin(l) = min(dTz(l,:)); %K

%calculate N2(z)
N2z = zeros(1,length(z));
N1z = zeros(1,length(z));
for r = 0:dr:core_radius-dr
    for i = 1:length(z)
        N2z(i) = N2z(i)+N2(Pp(i),Psforward(:,i).'+Psbackward(:,i).',r,1:length(lambdaS))/numPoints;
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
lambdaF_SL(l) = mean(lambdaFz_SL);


% dQz_max = h*(tauRad/tau*freqP-freqF_SE)*2./(cs_aP+cs_eP)*Am.*cs_aP*N0*core_radius^2/tauRad/w^2;
% part1 = h*(tauRad/tau*freqP-freqF_SE);
% dirPart1 = diff(part1(:))./diff(lambdaP(:));
% part2 = 2./(cs_aP+cs_eP)*Am.*cs_aP*N0*core_radius^2/tauRad/w^2;
% dirPart2 = diff(part2(:))./diff(lambdaP(:));
% 
% figure
% hold on
% plot(lambdaP*1e9,[NaN dirPart1.']);
% plot(lambdaP*1e9,[NaN dirPart2.']);

% %plot Qtot for different lengths
% Qtotz = zeros(1,length(z)-1);
% for zee = 2:length(z)
%     Qtotz(zee-1) = trapz(z(1:zee),dQz_dt2(l,1:zee));
% end
% figure(13)
% plot(z(2:end),Qtotz*1e3)
% xlabel('length (m)');
% ylabel('Total heat extracted (mW)');


% %graphs
% Pump power
figure(1)
grid on
hold on
plot(z,Pp,'DisplayName',sprintf('%g nm',lambdaP(l)*1e9));
xlabel('fiber length (m)');
ylabel('Pp(z) (W)');
title('Change in pump Power along fiber');

%forward signal
figure(2)
hold on
grid on
for i = validWavelengths
    plot(z,Psforward(i,:),'DisplayName',sprintf('%g nm',lambdaS(i)*1e9));
end
xlabel('fiber length (m)');
ylabel('Psforward(z) (W)');
title('Change in forward signal along fiber');
legend('show');

%backward signal
figure(3)
hold on
grid on
for i = validWavelengths
    plot(z,Psbackward(i,:),'DisplayName',sprintf('%g nm',lambdaS(i)*1e9));
end
xlabel('fiber length (m)');
ylabel('Psbackward(z) (W)');
title('Change in backward signal along fiber');
legend('show')
% 
% %forward Signal Spectrum
% figure(4)
% hold on
% grid on
% PsforwardNorm = Psforward(:,length(z))*1e3/dlam/1e9;
% plot(lambdaS*1e9,PsforwardNorm,'DisplayName',sprintf('%g nm, mean = %g nm',lambdaP(l)*1e9,lambdaFF_ASE(l)*1e9));
% xlabel('wavelength (nm)');
% ylabel('power (mW/nm)');
% title('Normalized Spectrum of forward signal')
% 
% %backward Signal Spectrum
% figure(5)
% hold on
% grid on
% PsbackwardNorm = Psbackward(:,1)*1e3/dlam/1e9;
% plot(lambdaS*1e9,PsbackwardNorm,'DisplayName',sprintf('%g nm, mean = %g nm',lambdaP(l)*1e9,lambdaFB_ASE(l)*1e9));
% xlabel('wavelength (nm)');
% ylabel('power (mW/nm)');
% title('Normalized Spectrum of backward signal')


% %N2(z)
% figure(6)
% grid on
% hold on
% plot(z,N2z/N0,'DisplayName',sprintf('%g nm',lambdaP(l)*1e9));
% %plot(z,N0*ones(length(z)))
% xlabel('fiber length (m)');
% ylabel('N2(z)/N0');
% title('N2 along fiber');

%dQ/dt as a function of z
figure(7)
%grid on
box on
hold on
plot(z,dQz_dt2(l,:)*1e3,'DisplayName',sprintf('%g nm',lambdaP(l)*1e9));
xlabel('Position z along fiber (m)');
ylabel('Heat extracted per unit length (mW/m)');
title('dQ/dt along the fiber');
% 
%dT as a function of z
figure(8)
hold on
box on
plot(z,dTz(l,:),'DisplayName',sprintf('%g nm',lambdaP(l)*1e9));
xlabel('Position z along the fiber (m)');
ylabel('Change in temperature (K)');
title('change in Temperature along the fiber');

end

% figure(1) %pump
% legend('show');
% figure(4) %forward ASE spectrum
% legend('show');
% figure(5) %backward ASE spectrum
% legend('show');
% figure(6) %N2(z)
% legend('show');
% figure(7) %dQ/dt(z)
% legend('show')
% figure(8) %dTz
% legend('show');
% 
% %dTmin vs. pump wavelength
% figure(9)
% hold on
% plot(lambdaP*1e9,dTmin);
% xlabel('Pump Wavelength (nm)');
% ylabel('Change in Temperature (K)');
% title('Temperature Change vs. Pump Wavelength');
% grid on

% %dQ/dt(z) vs. Pump wavelength
% figure(10)
% hold on
% plot(lambdaP*1e9,avgdQz_dt2(1,:));
% xlabel('Pump Wavelength (nm)');
% ylabel('dQ/dt (W/m)');
% title('dQ/dt vs. Pump Wavelength');
% grid on

%total heat extracted and efficiency vs. pump wavelength
figure(11)
hold on
yyaxis left
plot(lambdaP*1e9,dQ_dt2*1e3);
ylabel('Total heat extracted (mW)');
yyaxis right
plot(lambdaP*1e9,dQ_dt2./Pp0*100)
ylabel('Cooling efficiency (%)');
xlabel('Pump wavelength (nm)');
title('Total heat extracted and cooling efficiency vs. pump wavelength');
box on
pbaspect([2.73 1 1]);

%Pp0 vs. wavelength
figure(12)
hold on
plot(lambdaP*1e9,Pp0*1e3);
plot(lambdaP*1e9,Popt*1e3);
plot(lambdaP*1e9,Pp0.*(lambdaP-lambdaF_SE)./lambdaP)
xlabel('wavelength (nm)');
ylabel('input pump (mW)');
title('pump vs. wavelength');

%total heat extracted and Pp0 vs. pump wavelength
figure(13)
hold on
yyaxis left
plot(lambdaP*1e9,dQ_dt2*1e3);
ylabel('Total heat extracted (mW)');
yyaxis right
plot(lambdaP*1e9,Pp0*1e3)
ylabel('Optimum pump power (mW)');
xlabel('Pump wavelength (nm)');
title('Total heat extracted and cooling efficiency vs. pump wavelength');
box on

% 
% %mean flouresent wavelength for ASE vs. pump wavelength
% figure(11)
% hold on
% %plot(lambdaP*1e9,lambdaFF_ASE);
% plot(lambdaP*1e9,lambdaFB_ASE*1e9);
% xlabel('Pump Wavelength (nm)');
% ylabel('Mean flouresent wavelength (nm)');
% title('Mean wavelength of backward ASE vs. Pump Wavelength');
% grid on

% %gain of select wavelengths vs. pump wavelength
% figure(12)
% hold on
% for i = indexGain 
%     %plot(lambdaP*1e9,avgGain(:,i),'DisplayName',sprintf('%g nm',lambdaS(i)*1e9));
%     plot(lambdaP*1e9,gainASE(:,i),'DisplayName',sprintf('%g nm',lambdaS(i)*1e9));
% end
% xlabel('Pump Wavelength (nm)');
% ylabel('Gain (m^{-1})');
% title('Gain vs. Pump Wavelength');
% grid on
% legend('show');

% %gain of select wavelengths vs. pump wavelength
% %dQ/dt(z) vs. Pump wavelength (uses the avg value of dQ/dt)
% figure(13)
% yaxis = plotyy(lambdaP*1e9,avgGain(:,indexGain),lambdaP*1e9,avgdQz_dt2(1,:)*1e3);
% %yaxis = plotyy(lambdaP*1e9,avgGain(:,indexGain),lambdaP*1e9,dQ*1e3);
% gainAxis = ylim(yaxis(1));
% gainTicks = [-2*gainAxis(end)/3 -1*gainAxis(end)/3 0 gainAxis(end)/3 2*gainAxis(end)/3 gainAxis(end)];
% set(yaxis(1),'YTick',round(gainTicks,2),'ylim',[gainAxis(1) gainAxis(end)],'Box','off','ycolor','black');
% dQAxis = ylim(yaxis(2));
% dQTicks = [-2*dQAxis(end)/3 -1*dQAxis(end)/3 0 dQAxis(end)/3 2*dQAxis(end)/3 dQAxis(end)];
% set(yaxis(2),'YTick',round(dQTicks,2),'ylim',[dQAxis(1) dQAxis(end)],'Box','off','ycolor','black');
% ylabel(yaxis(1),'Gain (m^{-1})');
% ylabel(yaxis(2),'dQ/dt (mW/m)');
% xlabel('Pump Wavelength (nm)');
% title('Gain and dQ/dt vs. Pump Wavelength');
% legend(sprintf('%g nm',lambdaGains(1)*1e9),sprintf('%g nm',lambdaGains(2)*1e9),sprintf('%g nm',lambdaGains(3)*1e9),...
%     sprintf('%g nm',lambdaGains(4)*1e9),'dQ/dt');
% grid on


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


