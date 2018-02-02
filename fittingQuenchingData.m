close all 

%concentration
x = [1.15 1.53 1.90 2.3 4.6 6.5]*1e25; %1/m^3
%total lifteime
y = [780 800 820 840 805 615]*1e-3; %s

%define fit function
myfittype = fittype('c.*(1+a.*x)./(1+9/2/pi.*(x./b).^2)',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a','b','c'})

%do fit
myfit = fit(x',y',myfittype,'start',[11e-17 10e25 .75])

figure(1)
semilogx(x,y)
hold on
plot(myfit)

