%numerical result for change in power across a MMF for one pump and
%broadboand signal. Automatically converges. Does NOT approximate intensity in
%fiber as P(z)/Area.
%Automatically discards signal wavelengths for which gain is too low
%includes non-radiative term and loss

close all;
clear all;

len = .29;
dz = .001;
Pp0 = 71e-3;%[1:10:300];
%Psforward0 = 0;
%Psbackward0 = ones(21,1)*Pp0/21; %%GUESS THIS!!!!
maxIterations = 7;
maxTries = 0;
gainCutOff = 0.08; %if a signal wavelength has a gain factor that is less 
                   %than gainCutOff percent of the max gain factor, then
                   %the signal at that wavelength is set to 0;
minGain = -2500; %if a signal wavelength has a gain that is less 
               %than minGain than the signal at that wavelength is set to 0
maxError = 1e-8;

core_radius = 3.1e-6;%1.5e-6; 
Ac = pi*core_radius^2; %m^2, pump area
etta = 20; %fraction of the core diameter (radius?) at which the distribution is at full width
tauRad = 0.85e-3;
tauNonRad = .176e9;
tauRatio = tauRad/tauNonRad;
c = 3e8;
h = 6.63e-34; %J*s
sb = 5.67e-8; %W/m^2/K^4

%cross sectional areas 
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
lambda_f = trapz(cs_emsRAW(:,1),integral1)/trapz(cs_emsRAW(:,1),integral2)*1e-9;
freqF = c/lambda_f;


%create a variable to help me compare cross section values
crossSection = cs_abs;
crossSection(:,3) = cs_ems(:,2);
crossSection(:,4) = cs_abs(:,2)./cs_ems(:,2);

%waveleghth of laser/cooling pump
dlam = 10e-9; %resolution of the signal spectrum
lambdaP = 1020e-9;
lambdaS = [930e-9:dlam:1050e-9]; 
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

% ncore = 1.49;
% nclad = 1.4889;
% lam = 1000e-9;
% NA = .12; 
% V = 2*pi*core_radius/lam*NA;
% w = core_radius*(0.65+1.619/V^1.5+2.879/V^6);

N0 = 1.81e26; %m^-3
gamma = 1;
Am = pi*(etta*core_radius)^2;
loss_b = 0.01*log(10)/10; %m^-1
loss_bs = 1*loss_b/2;%m^-1
loss_ba = 1*loss_b/2;%0; %m^-1
IsatPa = h*freqP/cs_aP/tauRad;
PsatPa = IsatPa*Am;
IsatSa = h*freqS./cs_aS/tauRad;
IsatP = h*freqP/(cs_aP+cs_eP)/tauRad;
PsatP = IsatP*Am;
IsatS = h*freqS./(cs_aS+cs_eS)/tauRad;


fr = 1/Am;
%N2(z) - Pstot is 1xlength(lambdaS) array where each element is the sum of
%the backward and forward power for that wavelength at position 
N2 = @(Pp,Pstot,r,vw) (Pp*fr/IsatPa+sum(Pstot*fr./IsatSa(vw)))/(tauRatio+1+Pp*fr/IsatP+sum(Pstot*fr./IsatS(vw)))*N0; 
%signal gain coefficients
gammaS = @(Pp,Pstot,r,vw) N2(Pp,Pstot,r,vw)*(cs_eS(vw)+cs_aS(vw))-cs_aS(vw)*N0; %returns and array of size 1xlength(validWavelengths) 
gammaSE = @(Pp,Pstot,r,vw) cs_eS(vw)*N2(Pp,Pstot,r,vw); %returns and array of size 1xlength(validWavelengths)

comTerm1 = @(pump,Psf,Psb,vw) (pump/PsatPa+sum((Psf+Psb)./IsatSa(vw))/Am);
comTerm2 = @(pump,Psf,Psb,vw) (pump/PsatP+sum((Psf+Psb)./IsatS(vw))/Am);
%pump power - Psf and Psb are 1xlength(validWavelengths) arrays where each element
%is the forward and backward power (respectively) for that wavelegth and
%postion.
dPp_dz = @(z,pump,Psf,Psb,vw) -pump*N0/etta^2*((tauRatio+1)*cs_aP+cs_aP*sum((Psf+Psb)./IsatS(vw))/Am...
    -(cs_aP+cs_eP)*sum((Psf+Psb)./IsatSa(vw))/Am)/(tauRatio+1+pump/IsatP/Am+sum((Psf+Psb)./IsatS(vw))/Am)...
    -loss_b*pump/etta^2;
%forward signal power
dPsforward_dz = @(z,pump,Psf,Psb,vw) Psf*N0/etta^2.*(comTerm1(pump,Psf,Psb,vw)/(tauRatio+1+comTerm2(pump,Psf,Psb,vw))*(cs_eS(vw)+cs_eS(vw))...
    -cs_aS(vw))+2*h*freqS(vw).*dvS(vw).*cs_eS(vw)*N0/etta^2*comTerm1(pump,Psf,Psb,vw)/(tauRatio+1+comTerm2(pump,Psf,Psb,vw))...
    -loss_b/etta^2*Psf;
%backward signal power
dPsbackward_dz = @(z,pump,Psf,Psb,vw) -Psb*N0/etta^2.*(comTerm1(pump,Psf,Psb,vw)/(tauRatio+1+comTerm2(pump,Psf,Psb,vw))*(cs_eS(vw)+cs_eS(vw))...
    -cs_aS(vw))-2*h*freqS(vw).*dvS(vw).*cs_eS(vw)*N0/etta^2*comTerm1(pump,Psf,Psb,vw)/(tauRatio+1+comTerm2(pump,Psf,Psb,vw))...
    +loss_b/etta^2*Psb;


z = 0:dz:len;  
Pp = zeros(1,length(z)); 
Psforward = zeros(length(lambdaS),length(z)); 
Psbackward = zeros(length(lambdaS),length(z)); 

Ptotsignal = zeros(1,length(Pp0)); %stores total signal power emitted from fiber for each Pp0
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
dTz = zeros(length(Pp0),length(z));  %stores the change in temp at a location z along the fiber
dTzVac = zeros(length(Pp0),length(z));  %stores the change in temp at a location z along the fiber in a vaccuum
dTavg = zeros(1,length(Pp0)); %stores the average change in temperature across the fiber 
dTVacAvg = zeros(1,length(Pp0)); %stores the average change in temperature across the fiber in a vaccuum. 

for p = 1:length(Pp0)
    
validWavelengths = find(ones(1,length(lambdaS)));
    
Psbackward(:,1) = ones(length(lambdaS),1)*Pp0(p)/length(lambdaS)/2; %%GUESS THIS!!!!
Pp(1) = Pp0(p); %W  

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
awl = 1:length(lambdaS); %indexes of all wavelengths
for zz = 1:length(z)
    dP_dz(p,zz) = N0/etta^2*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl))...
        *Pp(zz)*(cs_aP+cs_eP)-N0/etta^2*Pp(zz)*cs_aP;
    
    dSF_dz(p,zz) = N0/etta^2*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl))...
        *(T2(Psforward(:,zz).')+2*h*sum(cs_eS(validWavelengths).*freqS(validWavelengths).*dvS(validWavelengths)))...
        -N0/etta^2*T1(Psforward(:,zz).');

    dSB_dz(p,zz) = -N0/etta^2*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)/(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl))...
        *(T2(Psbackward(:,zz).')+2*h*sum(cs_eS(validWavelengths).*freqS(validWavelengths).*dvS(validWavelengths)))...
        +N0/etta^2*T1(Psbackward(:,zz).');
    
    dQz_dtSE(p,zz) = Ac*N0*comTerm1(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl)...
        /(tauRatio+1+comTerm2(Pp(zz),Psforward(:,zz).',Psbackward(:,zz).',awl))*h*freqF/tauRad;
      
    dQz_dtPhononLoss(p,zz) = (Pp(zz)+sum(Psforward(:,zz)+Psbackward(:,zz)))*loss_ba/etta^2;
    
    dQz_dtScat(p,zz) = (Pp(zz)+sum(Psforward(:,zz)+Psbackward(:,zz)))*loss_bs/etta^2;
    
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
Ptotsignal(p) = sum(Psforward(:,length(z))-Psforward(:,1)+Psbackward(:,1)-Psbackward(:,length(z)));
Pabs(p) = Pp(1)-Pp(length(z));
dQ_dt1(p) = Pabs(p)-Ptotsignal(p)+dQ_dtSE(p)+dQ_dtScat(p);

%calculate change in temperature
b = 63.5e-6; %m
SA = 2*pi*b*len;
density = 2.65e3; %kg/m^3
cv = 680; %J/kg/K
mcv = pi*b^2*len*density*cv;
rTemp = 300;
emis_f = 1;
emis_c = 0.1;
radius_c = 2e-3; %m
chi = (1-emis_c)*emis_f*b/emis_c/radius_c;
thermalC = 81.4; %W/m^2/K
dTz(p,:) = dQz_dt2(p,:)/2/pi/b/thermalC; %K
dTzVac(p,:) = (rTemp^4+dQz_dt2(p,:)*(1+chi)/(sb*2*pi*b*emis_f)).^0.25;
dTVacAvg(p) = mean(dTzVac(p,:));
dTavg(p) = mean(dTz(p,:)); %K


%calculate N2(z)
N2z = zeros(1,length(z));
for r = 0:dr:core_radius-dr
    for i = 1:length(z)
        N2z(i) = N2z(i)+N2(Pp(i),Psforward(:,i).'+Psbackward(:,i).',r,1:length(lambdaS))/numPoints;
    end
end


%graphs
% Pump power
%approx = Pp0(p)*exp(-cs_aP*N0*z*(1-exp(-2*core_radius^2/w^2)));
figure(1)
grid on
hold on
plot(z,Pp);
%plot(z,log(approx));
xlabel('z (m)');
ylabel('Pp(z) (W)');
title('Change in pump Power along fiber');
legend('non approx','approx');

% %forward signal
% figure(2)
% hold on
% grid on
% for i = validWavelengths
%     plot(z,Psforward(i,:),'DisplayName',sprintf('%g nm',lambdaS(i)*1e9));
% end
% xlabel('z (m)');
% ylabel('Psforward(z) (W)');
% title('Change in forward signal along fiber');
% legend('show');
% 
% %backward signal
% figure(3)
% hold on
% grid on
% for i = validWavelengths
%     plot(z,Psbackward(i,:),'DisplayName',sprintf('%g nm',lambdaS(i)*1e9));
% end
% xlabel('z (m)');
% ylabel('Psbackward(z) (W)');
% title('Change in backward signal along fiber');
% legend('show')
% 
% %pump and signals
% figure(4)
% hold on 
% grid on
% plot(z,Pp);
% if length(lambdaS) > 1
%     plot(z,sum(Psforward));
%     plot(z,sum(Psbackward));
% else
%     plot(z,Psforward);
%     plot(z,Psbackward);
% end
% xlabel('z (m)');
% ylabel('Power (W)');
% title('Change in Power along fiber');
% legend('pump','forward signal','backward signal');
% 
% %forward Signal Spectrum
% figure(5)
% hold on
% grid on
% PsforwardNorm = Psforward(:,length(z))*1e3/dlam/1e9;
% plot(lambdaS*1e9,PsforwardNorm,'DisplayName',sprintf('%g mW, mean = %g nm',Pp0(p)*1e3,mlambda_f*1e9));
% xlabel('wavelength (nm)');
% ylabel('power (mW/nm)');
% title('Normalized Spectrum of forward signal')

% %backward Signal Spectrum
% figure(6)
% hold on
% grid on
% PsbackwardNorm = Psbackward(:,1)*1e3/dlam/1e9;
% plot(lambdaS*1e9,PsbackwardNorm,'DisplayName',sprintf('%g mW, mean = %g nm',Pp0(p)*1e3,mlambda_b*1e9));
% xlabel('wavelength (nm)');
% ylabel('power (mW/nm)');
% title('Normalized Spectrum of backward signal')

% %Signal Gain
% figure(7) 
% hold on
% grid on
% for i = 10%1:length(lambdaS)
%     plot(z,gain(i,:),'DisplayName',sprintf('%g nm',lambdaS(i)*1e9));
% end
% xlabel('z (m)');
% ylabel('gain(z) (1/m)');
% title('Signal gain along fiber');
% legend('show')

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

% %N2(z)
% figure(8)
% grid on
% hold on
% plot(z,N2z/N0,'DisplayName',sprintf('%g mW',Pp0(p)*1e3));
% %plot(z,N0*ones(length(z)))
% xlabel('z (m)');
% ylabel('N2(z)/N0');
% title('N2 along fiber');

% %dQ/dt as a function of z
% figure(9)
% grid on 
% hold on
% plot(z,dQz_dt2(p,:),'DisplayName',sprintf('%g mW',Pp0(p)*1e3));
% xlabel('z (m)');
% ylabel('dQ/dt (J/s)');
% title('dQ/dt along the fiber');

% %dT in a vacuum as a function of z 
% figure(14)
% hold on
% grid on
% plot(z,dTzVac(p,:)-rTemp,'DisplayName',sprintf('%g mW',Pp0(p)*1e3));
% xlabel('z (m)');
% ylabel('dT (K)');
% title('change in Temperature along the fiber in a vaccuum');

%dT in air as a function of z 
figure(15)
hold on
grid on
plot(z,dTz(p,:),'DisplayName',sprintf('%g mW',Pp0(p)*1e3));
xlabel('z (m)');
ylabel('dT (K)');
title('change in Temperature along the fiber in air');
end
% 
% figure(5) %forward signal spectrum
% legend('show')
% figure(6) %backward signal spectrum
% legend('show')
% figure(8) %N2(z)
% legend('show');
% figure(9) %dQ/dt(z)
% legend('show')
% plot(z(maxdQz_dt2(2,:)),maxdQz_dt2(1,:));

% %total signal power vs. absorped pump power
% figure(10)
% grid on
% plot(Pabs*1e3,(Ptotsignal-dQ_dtSE)*1e3)
% xlabel('Pabs (mW)');
% ylabel('total signal and SE power emmitted (mW)');
% title('total signal+SE power vs. absorbed pump power');

%Q
% figure(11)
% grid on
% hold on
% %plot(Pp0,dQ_dt1);
% plot(Pp0,dQ_dt2);
% xlabel('Pp(0) (W)');
% ylabel('dQ/dt (J/s)')
% title('Q vs Initial pump power');
% %legend('Pabs-Ps','differential equation');
% % 

% %average T in Vac vs. Pp0
% figure(16)
% hold on
% plot(Pp0,dTVacAvg-rTemp)
% xlabel('Pp0 (W)')
% ylabel('Average Change in Temp (K)')
% title('Average delT in vacuum vs. Pp0')
% grid on


% %SE and ASE contribution to Q
% figure(12)
% grid on
% hold on
% plot(Pp0,Ptotsignal);
% plot(Pp0,dQ_dtSE);
% xlabel('Pp(0) (W)');
% ylabel('dQ_SE/dt (J/s)')
% title('SE contribution to Q vs. Initial pump power');
% legend('ASE','other SE');





% me = (Pp(length(z))-Pp(1))/len(1);
% ma = (approx(length(z))-approx(1))/len(1);
% mratio = me/ma
% logFactor = (1+exp(-2*core_radius^2/w^2))/2
% slopefactor = log(logFactor)

tempurature = dTz(1,1)
power = Pp(end)
