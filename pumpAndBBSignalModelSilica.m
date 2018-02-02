%numerical result for change in power across fiber for one pump and
%broadboand signal. Automatically converges. Does NOT approximate intensity in
%fiber as P(z)/Area.
%Automatically discards signal wavelengths for which gain is too low
%includes non-radiative term and loss

%close all;%
clear all;

len = .0065;
dz = len/300;
Pp0 = 0.3e-3;%2332e-3;%.001:0.01:.6;%[.2 .175 .15 .124];%[0.001:0.01:0.7];%.05:0.05:1;%[.09:.03:2 3];
%Psforward0 = 0;
%Psbackward0 = ones(21,1)*Pp0/21; %%GUESS THIS!!!!
maxIterations = 25;
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
tauRad = 0.85e-3;
tauNonRad = 5e-3;%211.65e-3;%1e10;%
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
lambdaF_SE2 = lambdaF_SE/.996;
freqF_SE = c/lambdaF_SE;

%create a variable to help me compare cross section values
crossSection = cs_abs;
crossSection(:,3) = cs_ems(:,2);
crossSection(:,4) = cs_abs(:,2)./cs_ems(:,2);

%waveleghth of laser/cooling pump
dlam = 1e-9; %resolution of the signal spectrum
lambdaP = 1020e-9;
lambdaS = lambdaP;%[850e-9:dlam:1150e-9]; 
%validWavelengths = find(ones(1,length(lambdaS)));
freqP = c/lambdaP;
freqS = c./lambdaS;
dvS = zeros(1,length(lambdaS));%c./(lambdaS-dlam/2)-c./(lambdaS+dlam/2); %zeros(1,length(lambdaS));%
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


lam = 1000e-9;
NA = .13; 
%V = 2*pi*core_radius/lam*NA;
V = 2.44;
w = core_radius*(0.65+1.619/V^1.5+2.879/V^6);
etta = 1-exp(-2*core_radius^2/w^2);

N0_wtPercent = 2.71;%[0.00005:.00005:.001 .001:.005:.1 .1:.1:3 3:15];
N0 = N0_wtPercent/2.71*1.81e26;  %m^-3
f = 1/Area;
gamma = 1;
totLoss = 0.02; %dB/m
loss_b = totLoss/etta*log(10)/10; %m^-1
loss_bs = loss_b/4;%m^-1
loss_ba = 3*loss_b/4+43.8;%0; %m^-1
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
QmaxNoLoss2 = -c*h*(1/lambdaP-1/lambdaF_SE2)/(cs_aP+cs_eP)*Am*cs_aP*N0/tau*(-2*core_radius^2/w^2);

Popt = IsatP*pi/2*(-2*core_radius^2./log(1-etta)).*(-1*(2-etta)+sqrt(etta.^2+4*(1-etta)*cs_aP*N0/loss_ba*(tau/tauRad*lambdaP/lambdaF_SE-1)))/2./(1-etta);

%Pp0 = Popt;

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
dPp_dz = @(z,pump,Psf,Psb,vw) pump*cs_aP*N0/comTerm2(pump,Psf,Psb,vw)...
    *(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb,vw)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb,vw)))...
    +pump*N0*(cs_aP*sum((Psf+Psb)./IsatS(vw))/Am-(cs_aP+cs_eP)*sum((Psf+Psb)./IsatSa(vw))/Am)/comTerm2(pump,Psf,Psb,vw)...
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
    
validWavelengths = find(ones(1,length(lambdaS)));
    
% if sum(dvS)>0 %initialize the backward signal ONLY IF consider ASE in simulation
%     Psbackward(:,1) = ones(length(lambdaS),1)*Pp0(p)/length(lambdaS)/2; %%GUESS THIS!!!!
% end
%Psbackward(:,1) = Psbackward980(:,1)*percent;
Pp(1) = Pp0(p); %W  

iterations = 0;
tries = 0;


while iterations <= maxIterations
    
    %forward propagation
    dz = abs(dz);
    for i=1:(length(z)-1) % calculation loop
        
        if i == 384
            test = 0;
        end

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
        if Pp(end) > Pp(1)
            Pp(end) = Pp0(p)*exp(-cs_aP*N0*len*(1-exp(-2*core_radius^2/w^2)));
        end
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
    if (sum(abs(Psforward(validWavelengths,1)) > maxError) == 0 && abs(Pp(1)-Pp0(p)) < maxError) %%|| iterations > maxIterations
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
    Pp(1) = Pp0(p);
end

%change Psforward/Psbackward so that data for invalid wavelengths are set to zero.
%data for valid wavelengths is kept the same.
invalidWavelengths = ~ismember(1:length(lambdaS),validWavelengths);
Psforward(invalidWavelengths,:) = zeros(sum(invalidWavelengths),length(z));
Psbackward(invalidWavelengths,:) = zeros(sum(invalidWavelengths),length(z));

%calculate gain 
gain = zeros(length(lambdaS),length(z));
avgGain = zeros(1,length(lambdaS)); %stores the average gain across the fiber for each wavelength
maxGain = zeros(1,length(lambdaS)); %stores the maximum gain across the fiber for each wavelength
numPoints = 5;
dr = core_radius/numPoints;
for r = 0:dr:core_radius-dr %takes the average gain across radius of fiber
    for i = 1:length(z)
        gain(:,i) = gain(:,i)+gammaS(Pp(i),Psforward(:,i).'+Psbackward(:,i).',r,1:length(lambdaS)).'/numPoints;
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
%T1 = @(Ps) sum(Ps./IsatSa)*h*freqF_SE/tauRad;%sum(Ps.*cs_aS);%
%T2 = @(Ps) sum(Ps./IsatS)*h*freqF_SE/tauRad;%sum(Ps.*(cs_aS+cs_eS));%
T1 = @(Ps) sum(Ps.*cs_aS);%
T2 = @(Ps) sum(Ps.*(cs_aS+cs_eS));%
E = sum(cs_eS.*freqS.*dvS);
%find dQ/dt as a function of z
awl = 1:length(lambdaS); %indexes of all wavelengths
for zz = 1:length(z)
    dP_dz(p,zz) = -N0*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)...
        *Pp(zz)*(cs_aP+cs_eP)*((exp(-2*core_radius^2/w^2)-1)...
        -(1+tauRatio)*log((tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)))...
        /comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl))...
        -N0*Pp(zz)*cs_aP*(1-exp(-2*core_radius^2/w^2));
    
    dSF_dz(p,zz) = -N0*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)...
        *(T2(Psforward(:,zz).')+2*h*sum(cs_eS(validWavelengths).*freqS(validWavelengths).*dvS(validWavelengths)))*((exp(-2*core_radius^2/w^2)-1)...
        -(1+tauRatio)*log((tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)))...
        /comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl))...
        -N0*T1(Psforward(:,zz).')*(1-exp(-2*core_radius^2/w^2));

    dSB_dz(p,zz) = N0*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)...
        *(T2(Psbackward(:,zz).')+2*h*sum(cs_eS(validWavelengths).*freqS(validWavelengths).*dvS(validWavelengths)))*((exp(-2*core_radius^2/w^2)-1)...
        -(1+tauRatio)*log((tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)))...
        /comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl))...
        +N0*T1(Psbackward(:,zz).')*(1-exp(-2*core_radius^2/w^2));

    dQz_dtSE(p,zz) = -Am*N0*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)*h*freqF_SE/tauRad...
        *log((tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)));

    dQzcorr_dtSE(p,zz) = -4*h*N0*E*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)...
        *((exp(-2*core_radius^2/w^2)-1)...
        -(1+tauRatio)*log((tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)))...
        /comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl));
    
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

%transient response of temperature
dT_dtVac = @(t,T) (-sb*SA*(T^4-rTemp^4)/(1+chi)+dQz_dt2(p,1)*len)/mcv; %K/s
tempInVac_t = ode45(dT_dtVac,0:.01:100,rTemp);


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
lambdaF_SL(p) = mean(lambdaFz_SL);


%calculate propogation of other modes
% ettaOM = 0.007; %amount of power that is in the core from other modes
% K = 0.01; %percent of power from pump that goes into other modes. 
% loss = loss_b;
% dPom_dz = @(zi,pump) -cs_aP*N1z(zi)*pump*ettaOM+cs_eP*N2z(zi)*pump*ettaOM-loss*pump;
% 
% dz = abs(dz);
% Pom = K*Pp0*ones(1,length(z)); %initial value for the power in the other modes. 
% for i=1:(length(z)-1) % calculation loop
% 
%     RK1 = dPom_dz(i,Pom(i));
%     RK2 = dPom_dz(i+1,Pom(i)+0.5*dz*RK1);
%     RK3 = dPom_dz(i+1,Pom(i)+0.5*dz*RK2);
%     RK4 = dPom_dz(i+1,Pom(i)+0.5*dz*RK3);
% 
%     Pom(i+1) = Pom(i) + (1/6)*(RK1+2*RK2+2*RK3+RK4)*dz;
% end
% 
% plot(z,Pom)
% hold on
% plot(z(40:end),Pom(1)*exp(-cs_aP*N0*z(1:18)*ettaOM))
% xlabel('z (m)')
% ylabel('Power in other modes (W)')
% legend('Solution','Approximation')
% title('Power in other modes along fiber')
% 
% powerInCore = Pom(end)*ettaOM


% %graphs
%Pump power
approx = Pp0(p)*exp(-cs_aP*N0*z*(1-exp(-2*core_radius^2/w^2)));
figure(1)
grid on
hold on
plot(z,Pp);
%plot(z,log(approx));
xlabel('z (m)');
ylabel('Pp(z) (W)');
title('Change in pump Power along fiber');
legend('non approx','approx');
% 
%forward signal
figure(2)
hold on
grid on
for i = validWavelengths
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
for i = validWavelengths
    plot(z,Psbackward(i,:),'DisplayName',sprintf('%g nm',lambdaS(i)*1e9));
end
xlabel('z (m)');
ylabel('Backward ASE Power (W)');
title('Change in Backward ASE Signal Along The Fiber');
% legend('show')

%pump and signals
figure(4)
hold on 
grid on
%plot(z,Pp);
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

%forward Signal Spectrum
figure(5)
hold on
grid on
PsforwardNorm = Psforward(:,length(z))*1e3/dlam/1e9;
plot(lambdaS*1e9,PsforwardNorm,'DisplayName',sprintf('%g mW, mean = %g nm',Pp0(p)*1e3,lambdaFF_ASE*1e9));
xlabel('wavelength (nm)');
ylabel('power (mW/nm)');
title('Normalized Spectrum of forward signal')

%backward Signal Spectrum
figure(6)
hold on
grid on
PsbackwardNorm = Psbackward(:,1)*1e3/dlam/1e9;
plot(lambdaS*1e9,PsbackwardNorm,'DisplayName',sprintf('%g mW, mean = %g nm',Pp0(p)*1e3,lambdaFB_ASE*1e9));
xlabel('wavelength (nm)');
ylabel('power (mW/nm)');
title('Normalized Spectrum of backward signal')

%Signal Gain
figure(7) 
hold on
grid on
for i =1:length(lambdaS)
    plot(z,gain(i,:),'DisplayName',sprintf('%g nm',lambdaS(i)*1e9));
end
xlabel('z (m)');
ylabel('gain(z) (1/m)');
title('Signal gain along fiber');
legend('show')
% 
% %Signal Gain vs. Wavelength
% figure(13)
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
%plot(z,N0*ones(length(z)))
xlabel('z (m)');
ylabel('N2(z)/N0');
title('N2 along fiber');
% 
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
% % 
%dT as a function of z
figure(14)
hold on
grid on
plot(z,dTz(p,:),'DisplayName',sprintf('P_{pump} = %g mW',Pp0(p)*1e3));
xlabel('z (m)');
ylabel('Change In Temperature (K)');
title('Change in Temperature Along the Fiber in Air');

% %dT in a vacuum as a function of z 
% figure(14)
% hold on
% grid on
% plot(z,dTzVac(p,:)-rTemp,'DisplayName',sprintf('P_{pump} = %g mW',Pp0(p)*1e3));
% xlabel('z (m)');
% ylabel('Change in Temperature (K)');
% title('Change in Temperature along the Fiber in a Vaccuum');

% %transient response of Temp in vacuum
% figure(15)
% hold on
% grid on
% plot(tempInVac_t.x,tempInVac_t.y,'DisplayName',sprintf('%g mW',Pp0(p)*1e3));
% xlabel('t (s)');
% ylabel('Temperature (K)');
% title('Transient Temperature along the fiber in a vaccuum');

end

% figure(5) %forward signal spectrum
% legend('show')
% figure(6) %backward signal spectrum
% legend('show')
% figure(8) %N2(z)
% legend('show');
% figure(13) %dTz
% legend('show');
% figure(14) %dTzVac
% legend('show');
% figure(9) %dQ/dt(z)
% legend('show')
% plot(z(maxdQz_dt2(2,:)),maxdQz_dt2(1,:));

% %average T in Vac vs. Pp0
% figure(16)
% hold on
% plot(Pp0,dTVacAvg)
% xlabel('Pp0 (W)')
% ylabel('Average Temp (K)')
% title('Average T in vacuum vs. Pp0')
% grid on
% 
% %max T in Vac vs. Pp0
% figure(17)
% hold on
% plot(Pp0,dTVacMax-rTemp)
% xlabel('Pump Power (W)')
% ylabel('Maximum temperature change (K)')
% title('Maximum temperature change in a vacuum vs. pump power')
% grid on



% %total signal power vs. absorped pump power
% figure(10)
% grid on
% plot(Pabs*1e3,(Ptotsignal-dQ_dtSE)*1e3)
% xlabel('Pabs (mW)');
% ylabel('total signal and SE power emmitted (mW)');
% title('total signal+SE power vs. absorbed pump power');
% % 
% %Q and cooling eff
% figure(11)
% hold on
% yyaxis left
% plot(Pp0*1e3,dQ_dt2*1e3);
% ylabel('Total heat extracted (mW)');
% yyaxis right
% plot(Pp0*1e3,dQ_dt2./Pp0*100)
% ylabel('Cooling efficiency (%)');
% xlabel('Pump power (mW)');
% title('Total heat extracted and cooling efficiency vs. pump power');
% box on

% %max dQ/dt(z) vs. Pp0
% figure(12)
% grid on
% hold on
% plot(Pp0,maxdQz_dt2(1,:));
% xlabel('Input pump power (W)');
% ylabel('Maximum heat extracted (J/s/m)')
% title('Maximum heat extracted in a vacuum vs input pump power');


% %SE and ASE contribution to Q
% figure(13)
% grid on
% hold on
% plot(Pp0,Ptotsignal);
% %plot(Pp0,dQ_dtSE);
% xlabel('Pp(0) (W)');
% ylabel('total Signal (W)')
% title('Total signal vs. Initial pump power');
% %legend('ASE','other SE');





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

maxLoss = cs_aP*N0*(tau/tauRad*freqF_SE/freqP-1); %1/m
maxCooling = -(tauRad/tau*h*freqP-h*freqF_SE)/(cs_eP+cs_aP)*Am*cs_aP*N0/tauRad*(-2*core_radius^2/w^2);

%find pump power such that the fiber has zero net temp change
findLossApprox = @(loss) maxCooling*0.9/PsatP+loss*etta+maxLoss*(log(maxLoss./loss)-1);
coolingApprox = @(loss) loss*etta+maxLoss*(log(maxLoss./loss)-1);
optPumpApprox = @(loss) PsatP*(-1+maxLoss./loss);

optPump = @(loss) IsatP*pi/2*(-2*core_radius^2/log(1-etta))*(-1*(2-etta)+sqrt(etta^2+4*(1-etta)*cs_aP*N0./loss*(tau/tauRad*lambdaP/lambdaF_SE-1)))/2/(1-etta);

cooling2 = @(loss) loss*etta.*optPump(loss)+PsatP*maxLoss*log((1+optPump(loss)/PsatP*(1-etta))./(1+optPump(loss)/PsatP));

cooling = @(loss) loss*etta*PsatP.*(sqrt(etta^2+4*(1-etta)*maxLoss./loss)-(2-etta))/2/(1-etta)+PsatP*maxLoss...
    *log((1+(sqrt(etta^2+4*(1-etta)*maxLoss./loss)-(2-etta))/2)./(1+(sqrt(etta^2+4*(1-etta)*maxLoss./loss)-(2-etta))/2/(1-etta)));

findLossExact = @(loss) loss*etta*PsatP.*(sqrt(etta^2+4*(1-etta)*maxLoss./loss)-(2-etta))/2/(1-etta)+PsatP*maxLoss...
    *log((1+(sqrt(etta^2+4*(1-etta)*maxLoss./loss)-(2-etta))/2)./(1+(sqrt(etta^2+4*(1-etta)*maxLoss./loss)-(2-etta))/2/(1-etta)))-maxCooling*0.9;


% approxLoss90 = fzero(findLossApprox,0.1*log(10)/10)*10/log(10)*etta
% if isnan(approxLoss90)
%     approxLoss90 = fzero(findLossApprox,0.1*log(10)/10)*10/log(10)*etta
% end
% 
% exactLoss90 = fzero(findLossExact,0.001*log(10)/10)*10/log(10)*etta %dB/m of measured loss
% if isnan(exactLoss90)
%     exactLoss90 = fzero(findLossExact,0.0007)*10/log(10)*etta
% end


