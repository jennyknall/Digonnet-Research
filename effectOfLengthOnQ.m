%numerical result for change in power across fiber for one pump and
%broadboand signal. Automatically converges. Does NOT approximate intensity in
%fiber as P(z)/Area.
%Automatically discards signal wavelengths for which gain is too low
%includes non-radiative term and loss

close all;
clear all;


len = .3;%[0.1 .2 .3];%.1:.1:1;
numdz = 100;
%dz = .0001;
Pp0 = .5;%.16;%[.09:.03:1 2 3];
%Psforward0 = 0;
%Psbackward0 = ones(21,1)*Pp0/21; %%GUESS THIS!!!!
maxIterations = 1;
maxTries = 0;
gainCutOff = 0.08; %if a signal wavelength has a gain factor that is less 
                   %than gainCutOff percent of the max gain factor, then
                   %the signal at that wavelength is set to 0;
minGain = -2500; %if a signal wavelength has a gain that is less 
               %than minGain than the signal at that wavelength is set to 0
maxError = 1e-5;

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

%waveleghth of laser/cooling pump
dlam = 3e-9; %resolution of the signal spectrum
lambdaP = 976e-9;
lambdaS = [850e-9:dlam:1150e-9]; 
%validWavelengths = find(ones(1,length(lambdaS)));
freqP = c/lambdaP;
freqS = c./lambdaS;
%dvS = zeros(1,length(lambdaS));
dvS = c./(lambdaS-dlam/2)-c./(lambdaS+dlam/2); %zeros(1,length(lambdaS));%
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


lam = 1020e-9;
NA = .13; 
V = 2*pi*core_radius/lam*NA;
w = core_radius*(0.65+1.619/V^1.5+2.879/V^6);


N0 = 1.81e26; %m^-3
f = 1/Area;
gamma = 1;
loss_b = .0023; %m^-1
loss_bs = loss_b/2;%m^-1
loss_ba = loss_b/2;%0; %m^-1
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
dPsforward_dz = @(z,pump,Psf,Psb,vw) -Psf.*(cs_eS(vw)+cs_aS(vw))*N0*comTerm1(pump,Psf,Psb,vw)/comTerm2(pump,Psf,Psb,vw)*((exp(-2*core_radius^2/w^2)-1)...
    -(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb,vw)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb,vw)))/comTerm2(pump,Psf,Psb,vw))...
    -Psf.*cs_aS(vw)*N0*(1-exp(-2*core_radius^2/w^2))-2*h*freqS(vw).*dvS(vw).*cs_eS(vw)*N0*comTerm1(pump,Psf,Psb,vw)/comTerm2(pump,Psf,Psb,vw)...
    *((exp(-2*core_radius^2/w^2)-1)-(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb,vw)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb,vw)))/comTerm2(pump,Psf,Psb,vw))...
    -loss_b*Psf*(1-exp(-2*core_radius^2/w^2));
%backward signal power
dPsbackward_dz = @(z,pump,Psf,Psb,vw) Psb.*(cs_eS(vw)+cs_aS(vw))*N0*comTerm1(pump,Psf,Psb,vw)/comTerm2(pump,Psf,Psb,vw)*((exp(-2*core_radius^2/w^2)-1)...
    -(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb,vw)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb,vw)))/comTerm2(pump,Psf,Psb,vw))...
    +Psb.*cs_aS(vw)*N0*(1-exp(-2*core_radius^2/w^2))+2*h*freqS(vw).*dvS(vw).*cs_eS(vw)*N0*comTerm1(pump,Psf,Psb,vw)/comTerm2(pump,Psf,Psb,vw)...
    *((exp(-2*core_radius^2/w^2)-1)-(1+tauRatio)*log((tauRatio+1+comTerm2(pump,Psf,Psb,vw)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(pump,Psf,Psb,vw)))/comTerm2(pump,Psf,Psb,vw))...
    +loss_b*Psb*(1-exp(-2*core_radius^2/w^2));



%z = 0:dz:len(l);  
Pp = zeros(1,numdz); 
Psforward = zeros(length(lambdaS),numdz); 
Psbackward = zeros(length(lambdaS),numdz); 

Ptotsignal = zeros(1,length(len)); %stores total signal power emitted from fiber for each len
Pabs = zeros(1,length(len)); %stores values for pump power absorbed for each len
dQ_dt1 = zeros(1,length(len)); %stores values for dQ/dt for each len
dQz_dt2 = zeros(length(len),numdz); %stores dQ/dt as a function of z for each len
dQz_dtSE = zeros(length(len),numdz); %stores the SE contribution to dQ/dt as a function of z for each len
dQzcorr_dtSE = zeros(length(len),numdz); %stores the correction term the SE cooling as a funciton of z for each Pp0
dQz_dtPhononLoss = zeros(length(len),numdz); %stores the phonon loss contribution to dQ/dt as a function of z for each len
dQz_dtScat = zeros(length(len),numdz); %stores the scattering loss contribution to dQ/dt as a function of z for each len
dSF_dz = zeros(length(len),numdz); %stores the forward signal contribution to dQ/dt as a function of z for each len
dSB_dz = zeros(length(len),numdz); %stores the backward signal contribution to dQ/dt as a function of z for each len 
dP_dz = zeros(length(len),numdz);%stores the Pump abs contribution to dQ/dt as a function of z for each len
dQ_dt2 = zeros(1,length(len)); %stores values for dQ/dt for each len
dQ_dtSE = zeros(1,length(len)); %stores the SE contribution to dQ/dt for each len
dQ_dtSEcorr = zeros(1,length(len)); %stores the SE correction term to dQ/dt for each Pp0
dQ_dtPhononLoss = zeros(1,length(len)); %stores the Phonon loss contribution to dQ/dt for each len
dQ_dtScat = zeros(1,length(len)); %stores the scattering loss contribution to dQ/dt for each len
dQ_dtSF = zeros(1,length(len)); %stores the forward signal contribution to dQ/dt for each len
dQ_dtSB = zeros(1,length(len)); %stores the backward signal contribution to dQ/dt for each len
dQ_dtP = zeros(1,length(len)); %stores the pump abs contribution to dQ/dt for each len
maxdQz_dt2 = zeros(2,length(len)); %stores the coordinates of the maximum value of dQz_dt2
                                    %first row: max values
                                    %second row: index of z where max value is
lambdaF_SL = zeros(1,length(len)); %stores the mean flouresnce wavelength for each pump power calculated from the 
                                    %populations of each sublevel.
dTz = zeros(length(len),numdz);  %stores the change in temp at a location z along the fiber
dTavg = zeros(1,length(len)); %stores the average change in temperature across the fiber 
dTmin = zeros(1,length(len)); %stores the minimum change in temperature across the fiber (most cooling if negative)

for l = 1:length(len)

dz = len(l)/numdz;
z = 0:dz:len(l)-dz;
    
validWavelengths = find(ones(1,length(lambdaS)));
    
Psbackward(:,1) = ones(length(lambdaS),1)*Pp0/length(lambdaS)/2; %%GUESS THIS!!!!
Pp(1) = Pp0; %W  

iterations = 0;
tries = 0;


while iterations <= maxIterations
 
    %forward propagation
    dz = abs(dz);
    for i=1:(length(z)-1) % calculation loop

        RK1p = dPp_dz(z(i),Pp(i),Psforward(validWavelengths,i).',Psbackward(validWavelengths,i).',validWavelengths);
        RK1f = dPsforward_dz(z(i),Pp(i),Psforward(validWavelengths,i).',Psbackward(validWavelengths,i).',validWavelengths); %array of 1xlength(lambdaS)
        RK1b = dPsbackward_dz(z(i),Pp(i),Psforward(validWavelengths,i).',Psbackward(validWavelengths,i).',validWavelengths); %array of 1xlength(lambdaS)

        RK2p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(validWavelengths,i).'+0.5*dz*RK1f,Psbackward(validWavelengths,i).'+0.5*dz*RK1b,validWavelengths);
        RK2f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(validWavelengths,i).'+0.5*dz*RK1f,Psbackward(validWavelengths,i).'+0.5*dz*RK1b,validWavelengths);
        RK2b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(validWavelengths,i).'+0.5*dz*RK1f,Psbackward(validWavelengths,i).'+0.5*dz*RK1b,validWavelengths);

        RK3p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(validWavelengths,i).'+0.5*dz*RK2f,Psbackward(validWavelengths,i).'+0.5*dz*RK2b,validWavelengths);
        RK3f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(validWavelengths,i).'+0.5*dz*RK2f,Psbackward(validWavelengths,i).'+0.5*dz*RK2b,validWavelengths);
        RK3b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(validWavelengths,i).'+0.5*dz*RK2f,Psbackward(validWavelengths,i).'+0.5*dz*RK2b,validWavelengths);

        RK4p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(validWavelengths,i).'+0.5*dz*RK3f,Psbackward(validWavelengths,i).'+0.5*dz*RK3b,validWavelengths);
        RK4f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(validWavelengths,i).'+0.5*dz*RK3f,Psbackward(validWavelengths,i).'+0.5*dz*RK3b,validWavelengths);
        RK4b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(validWavelengths,i).'+0.5*dz*RK3f,Psbackward(validWavelengths,i).'+0.5*dz*RK3b,validWavelengths);

        Pp(i+1) = Pp(i) + (1/6)*(RK1p+2*RK2p+2*RK3p+RK4p)*dz;
        Psforward(validWavelengths,i+1) = Psforward(validWavelengths,i) + (1/6)*(RK1f+2*RK2f+2*RK3f+RK4f).'*dz;  % main equation
        Psbackward(validWavelengths,i+1) = Psbackward(validWavelengths,i) + (1/6)*(RK1b+2*RK2b+2*RK3b+RK4b).'*dz;  % main equation
    end
    
    %check/modify boundary conditions
    if sum(abs(Psbackward(validWavelengths,length(z))) > maxError) == 0
        break
    else
        Psbackward(validWavelengths,length(z)) = 0;
    end
    
    %backward propagation
    dz = -abs(dz);
    for i=flip(2:length(z)) % calculation loop
        
        RK1p = dPp_dz(z(i),Pp(i),Psforward(validWavelengths,i).',Psbackward(validWavelengths,i).',validWavelengths);
        RK1f = dPsforward_dz(z(i),Pp(i),Psforward(validWavelengths,i).',Psbackward(validWavelengths,i).',validWavelengths); %array of 1xlength(lambdaS)
        RK1b = dPsbackward_dz(z(i),Pp(i),Psforward(validWavelengths,i).',Psbackward(validWavelengths,i).',validWavelengths); %array of 1xlength(lambdaS)

        RK2p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(validWavelengths,i).'+0.5*dz*RK1f,Psbackward(validWavelengths,i).'+0.5*dz*RK1b,validWavelengths);
        RK2f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(validWavelengths,i).'+0.5*dz*RK1f,Psbackward(validWavelengths,i).'+0.5*dz*RK1b,validWavelengths);
        RK2b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK1p,Psforward(validWavelengths,i).'+0.5*dz*RK1f,Psbackward(validWavelengths,i).'+0.5*dz*RK1b,validWavelengths);

        RK3p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(validWavelengths,i).'+0.5*dz*RK2f,Psbackward(validWavelengths,i).'+0.5*dz*RK2b,validWavelengths);
        RK3f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(validWavelengths,i).'+0.5*dz*RK2f,Psbackward(validWavelengths,i).'+0.5*dz*RK2b,validWavelengths);
        RK3b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK2p,Psforward(validWavelengths,i).'+0.5*dz*RK2f,Psbackward(validWavelengths,i).'+0.5*dz*RK2b,validWavelengths);

        RK4p = dPp_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(validWavelengths,i).'+0.5*dz*RK3f,Psbackward(validWavelengths,i).'+0.5*dz*RK3b,validWavelengths);
        RK4f = dPsforward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(validWavelengths,i).'+0.5*dz*RK3f,Psbackward(validWavelengths,i).'+0.5*dz*RK3b,validWavelengths);
        RK4b = dPsbackward_dz(z(i)+0.5*dz,Pp(i)+0.5*dz*RK3p,Psforward(validWavelengths,i).'+0.5*dz*RK3f,Psbackward(validWavelengths,i).'+0.5*dz*RK3b,validWavelengths);

        Pp(i-1) = Pp(i) + (1/6)*(RK1p+2*RK2p+2*RK3p+RK4p)*dz;
        Psforward(validWavelengths,i-1) = Psforward(validWavelengths,i) + (1/6)*(RK1f+2*RK2f+2*RK3f+RK4f).'*dz;  % main equation
        Psbackward(validWavelengths,i-1) = Psbackward(validWavelengths,i) + (1/6)*(RK1b+2*RK2b+2*RK3b+RK4b).'*dz;  % main equation
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
        gainS = zeros(length(lambdaS),length(z));
        gainSE = zeros(length(lambdaS),length(z));
        avgGainS = zeros(1,length(lambdaS)); 
        avgGainSE = zeros(1,length(lambdaS));
        for i = 1:length(z)
            gainS(:,i) = gainS(:,i)+gammaS(Pp(i),Psforward(:,i).'+Psbackward(:,i).',0,1:length(lambdaS)).';
            gainSE(:,i) = gainSE(:,i)+(gammaSE(Pp(i),Psforward(:,i).'+Psbackward(:,i).',0,1:length(lambdaS))*2*h.*freqS.*dvS*fr(0)).';
        end
        for i = 1:length(validWavelengths)
            avgGainS(i) = mean(gainS(i,:));
            avgGainSE(i) = mean(gainSE(i,:));
        end
        gainFactor = avgGainSE./avgGainS;
        newValidWavelengths = find(avgGainS > minGain);
        %newValidWavelengths = validWavelengths((abs(gainFactor) > gainCutOff*max(gainFactor)));
        if isequal(newValidWavelengths,validWavelengths)
            break
        else
            validWavelengths = newValidWavelengths;
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
T1 = @(Ps) sum(Ps.*cs_aS);
T2 = @(Ps) sum(Ps.*(cs_aS+cs_eS));
E = sum(cs_eS.*freqS.*dvS);
%find dQ/dt as a function of z
awl = 1:length(lambdaS); %indexes of all wavelengths
for zz = 1:length(z)
    dP_dz(l,zz) = -N0*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)...
        *Pp(zz)*(cs_aP+cs_eP)*((exp(-2*core_radius^2/w^2)-1)...
        -(1+tauRatio)*log((tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)))...
        /comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl))...
        -N0*Pp(zz)*cs_aP*(1-exp(-2*core_radius^2/w^2));
    
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
    dQzcorr_dtSE(l,zz) =0;
%     dQzcorr_dtSE(l,zz) = -4*h*N0*E*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)/comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)...
%         *((exp(-2*core_radius^2/w^2)-1)...
%         -(1+tauRatio)*log((tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)*exp(-2*core_radius^2/w^2))/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)))...
%         /comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl));
%     
    dQz_dtPhononLoss(l,zz) = (Pp(zz)+sum(Psforward(:,zz)+Psbackward(:,zz)))*loss_ba*(1-exp(-2*core_radius^2/w^2));
    
    dQz_dtScat(l,zz) = (Pp(zz)+sum(Psforward(:,zz)+Psbackward(:,zz)))*loss_bs*(1-exp(-2*core_radius^2/w^2));
    
    dQz_dt2(l,zz) = -dP_dz(l,zz)-dSF_dz(l,zz)+dSB_dz(l,zz)-dQz_dtSE(l,zz)+dQzcorr_dtSE(l,zz)+dQz_dtPhononLoss(l,zz);
end
[maxdQz_dt2(1,l), maxdQz_dt2(2,l)] = max(dQz_dt2(l,:));
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
b = 62.5e-6; %m
thermalC = 81.4; %W/m^2/K
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


%graphs
% % Pump power
% figure(1)
% grid on
% hold on
% plot(0:numdz-1,Pp,'DisplayName',sprintf('%g m',len(l)));
% xlabel('normalized fiber length');
% ylabel('Pp(z) (W)');
% title('Change in pump Power along fiber');

% %forward signal
figure(2)
hold on
grid on
for i = validWavelengths
    plot(0:numdz-1,Psforward(i,:),'DisplayName',sprintf('%g nm',lambdaS(i)*1e9));
end
xlabel('normalized fiber length');
ylabel('Psforward(z) (W)');
title('Change in forward signal along fiber');
legend('show');

%backward signal
figure(3)
hold on
grid on
for i = validWavelengths
    plot(0:numdz-1,Psbackward(i,:),'DisplayName',sprintf('%g nm',lambdaS(i)*1e9));
end
xlabel('normalized fiber length');
ylabel('Psbackward(z) (W)');
title('Change in backward signal along fiber');
legend('show')

% 
% %N2(z)
% figure(4)
% grid on
% hold on
% plot(0:numdz-1,N2z/N0,'DisplayName',sprintf('%g m',len(l)));
% %plot(z,N0*ones(length(z)))
% xlabel('normalized fiber length');
% ylabel('N2(z)/N0');
% title('N2 along fiber');

%dQ/dt as a function of z
figure(5)
grid on 
hold on
plot(0:numdz-1,dQz_dt2(l,:),'DisplayName',sprintf('%g m',len(l)));
xlabel('normalized fiber length');
ylabel('dQ/dt (J/s)');
title('dQ/dt along the fiber');

% %dT as a function of z
% figure(6)
% hold on
% grid on
% plot(0:numdz-1,dTz(l,:),'DisplayName',sprintf('%g m',len(l)));
% xlabel('Normalized Fiber Length (%)');
% ylabel('Change in temperature (K)');
% title('Change in Temperature Along The Fiber');

end

% figure(1) %pump
% legend('show');
% figure(4) %N2(z)
% legend('show');
% figure(5) %dQ/dt(z)
% legend('show')
% figure(6) %dTz
% legend('show');

% %Q
% figure(7)
% grid on
% hold on
% %plot(Pp0,dQ_dt1);
% plot(Pp0,dQ_dt2);
% xlabel('Pp(0) (W)');
% ylabel('dQ/dt (J/s)')
% title('Q vs Initial pump power');
% %legend('Pabs-Ps','differential equation');
% 
% %dTmin vs. length of fiber
% figure(8)
% hold on
% plot(len,dTmin);
% xlabel('Length of Fiber (m)');
% ylabel('Change in Temperature (K)');
% title('Temperature Change vs. Length of Fiber');
% grid on
% 
% %SE and ASE contribution to Q
% figure(9)
% grid on
% hold on
% plot(len,Ptotsignal);
% %plot(Pp0,dQ_dtSE);
% xlabel('Fiber Length (m)');
% ylabel('total Signal (W)')
% title('Total signal vs. fiber length');
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


