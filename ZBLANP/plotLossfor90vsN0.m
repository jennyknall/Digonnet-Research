%close all
clear all

N0_wtPercent = [0.00005:.00005:.001 .001:.005:.1 .1:.01:4];
conc = N0_wtPercent*2.42e26;  %m^-3

percentCooling = 0.9;

core_radius = 3.1e-6; 
tauRad = 1.7e-3;
Nc = 6.47e27; %m^-3
c = 3e8;
h = 6.63e-34; %J*s
kT = 4.11e-21; %J
sb = 5.67e-8; %W/m^2/K^4

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

%waveleghth of laser/cooling pump
lambdaP = 1015e-9;
freqP = c/lambdaP;
%cross sectional areas for a given wavelength
indexP = find(round(cs_abs(:,1)*1e9 - lambdaP*1e9) == 0);
cs_aP = cs_abs(indexP,2);
cs_eP = cs_ems(indexP,2);


lam = 1000e-9;
NA = .13; 
V = 2*pi*core_radius/lam*NA;
w = core_radius*(0.65+1.619/V^1.5+2.879/V^6);
etta = 1-exp(-2*core_radius^2/w^2);
Am = pi*w^2/2;

IsatPa = h*freqP/cs_aP/tauRad;
PsatPa = IsatPa*pi*w^2/2;
IsatP = h*freqP/(cs_aP+cs_eP)/tauRad;
PsatP = IsatP*pi*w^2/2;

loss90 = zeros(1,length(conc));
loss90approx = zeros(1,length(conc));
lossApproxMinus = zeros(1,length(conc));
lossApproxPlus = zeros(1,length(conc));

approxFactor = zeros(1,length(conc));

for i = 1:length(conc)
    
N0 = conc(i);

tauNonRad = 2*pi/9*tauRad*(Nc/N0)^2;%67.6e-3;%211.65e-3;%1e10;%
tauRatio = tauRad/tauNonRad;
tau = tauRad*tauNonRad/(tauRad+tauNonRad);

maxLoss = cs_aP*N0*(tau/tauRad*freqF_SE/freqP-1); %1/m
maxCooling = -(tauRad/tau*h*freqP-h*freqF_SE)/(cs_eP+cs_aP)*Am*cs_aP*N0/tauRad*(-2*core_radius^2/w^2);

%find pump power such that the fiber has zero net temp change

optPump = @(loss) IsatP*pi/2*(-2*core_radius^2/log(1-etta))*(-1*(2-etta)+sqrt(etta^2+4*(1-etta)*cs_aP*N0./loss*(tau/tauRad*lambdaP/lambdaF_SE-1)))/2/(1-etta);

%uses optPump function
cooling2 = @(loss) loss*etta.*optPump(loss)+PsatP*maxLoss*log((1+optPump(loss)/PsatP*(1-etta))./(1+optPump(loss)/PsatP));

%doesn't use optPump function - gives same result as cooling2
cooling = @(loss) loss*etta*PsatP.*(sqrt(etta^2+4*(1-etta)*maxLoss./loss)-(2-etta))/2/(1-etta)+PsatP*maxLoss...
    *log((1+(sqrt(etta^2+4*(1-etta)*maxLoss./loss)-(2-etta))/2)./(1+(sqrt(etta^2+4*(1-etta)*maxLoss./loss)-(2-etta))/2/(1-etta)));

findLossExact = @(loss) loss*etta*PsatP.*(sqrt(etta^2+4*(1-etta)*maxLoss./loss)-(2-etta))/2/(1-etta)+PsatP*maxLoss...
    *log((1+(sqrt(etta^2+4*(1-etta)*maxLoss./loss)-(2-etta))/2)./(1+(sqrt(etta^2+4*(1-etta)*maxLoss./loss)-(2-etta))/2/(1-etta)))-maxCooling*percentCooling;

findLossApprox = @(loss) (percentCooling*2*core_radius^2/w^2+0.948*log(1-etta))*2*(1-etta)/etta*maxLoss-(2-etta)*loss...
    +sqrt(etta^2*loss.^2+4*(1-etta)*loss*maxLoss);

D = percentCooling*2*core_radius^2/w^2+0.948*log(1-etta);
Kplus = 1/2+(2-etta)/2/etta*D+1/2*sqrt(1+2*(2-etta)/etta*D+D^2);
Kminus(i) = 1/2+(2-etta)/2/etta*D-1/2*sqrt(1+2*(2-etta)/etta*D+D^2);
lossApproxPlus(i) = Kplus*maxLoss*10/log(10)*etta;
lossApproxMinus(i) = Kminus(i)*maxLoss*10/log(10)*etta;

guess = maxLoss/100;

loss90(i) = fzero(findLossExact,guess)*10/log(10)*etta; %dB/m of measured loss
if isnan(loss90(i))
    loss90(i) = fzero(findLossExact,guess/10)*10/log(10)*etta;
end

loss90approx(i) = fzero(findLossApprox,guess)*10/log(10)*etta; %dB/m of measured loss
if isnan(loss90approx(i))
    loss90approx(i) = fzero(findLossApprox,guess/10)*10/log(10)*etta;
end

loss = loss90(i)/etta*log(10)/10;
pump = optPump(loss);
x = -pump/PsatP./(1+pump/PsatP)*etta;

degree = 50;
for n = 1:degree
    approxFactor(i) = approxFactor(i) + (-1)^(n+1)*x^(n)/n/log(1-etta);
end

end


figure(2)
hold on
box on
plot(N0_wtPercent,loss90*1e3)
%plot(N0_wtPercent,lossApproxMinus)
%plot(N0_wtPercent,(lossApproxMinus-loss90)./loss90)
%plot(N0_wtPercent,loss90approx)
%plot(N0_wtPercent,(loss90approx-loss90)./loss90)
xlabel('Yb concentration (wt% Yb)');
ylabel('Absorptive loss for 90% cooling (dB/km)')
%legend('exact','approx','error')
%ylim([0 max(loss90)])

figure(3)
box on
plot(N0_wtPercent,approxFactor)
xlabel('Yb concentration (wt% Yb_2O_3)');
ylabel('Approximation Factor \epsilon')


% loss = loss90(1)/etta*log(10)/10;%maxLoss/1000:maxLoss/1000:maxLoss;
% pump = optPump(loss);
% x = -pump/PsatP./(1+pump/PsatP)*etta;
% error = (x-x^2/2+x^3/3-x^4/4+x^5/5-x^6/6+x^7/7-x^8/8+x^9/9-x^10/10+x^11/11-x^12/12+x^13/13-x^14/14+...
%     x^15/15-x^16/16+x^17/17-x^18/18+x^19/19-x^20/20+x^21/21-x^22/22)/log(1-etta)
% 
% degree = 1000;
% error = 0;
% for n = 1:degree
%     error = error + (-1)^(n+1)*x^(n)/n/log(1-etta);
% end
% error
% 
% error = zeros(1,length(degree));
% for n = 2:degree
%     error(n) = error(n-1)+(-1)^(n)*x^(n-1)/(n-1)/log(1-etta);
% end


% findLossApprox = @(loss) (1.8*core_radius^2/w^2+error(end)*log(1-etta))*2*(1-etta)/etta*maxLoss-(2-etta)*loss...
%     +sqrt(etta^2*loss.^2+4*(1-etta)*loss*maxLoss);


