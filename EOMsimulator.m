clear all;
close all;

c = 3e8;

beta = 20;
fMod = .5;
f0 = 30;
A = 1;
B = 1;

dt = 1/f0/50;
t = 0:dt:1/fMod;
y1 = A*cos(2*pi*f0*t+beta*sin(2*pi*fMod*t));
% figure(1)
% hold on
% plot(t,y1);



% dtPerPeriod = 100;
% dt = 1/f0/dtPerPeriod;
Fs = 1/dt;
periodsPerN = 1;
% samples = periodsPerN*dtPerPeriod;
samples = 50*periodsPerN;

figure(1)
hold on

figure(2)
hold on

figure(3)
hold on

nTime = 200;
fmax = zeros(1,nTime);

for n = 1:nTime
    t = n*dt*(samples-1):dt:(n+1)*dt*(samples-1);
    y = A*cos(2*pi*f0*t+beta*sin(2*pi*fMod*t));
    figure(1)
    plot(y);
    S = fft(y);
    P2 = abs(S/samples);
    P1 = P2(1:samples/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(samples/2))/samples;
    figure(2)
    plot((f-f0)/fMod,P1)
    figure(3)
    plot(f,P1)
    [x,findex] = max(P1);
    fmax(n) = f(findex);
end

figure(4)
hold on
plot(fmax)
ylabel('Modulated frequency (Hz)')
xlabel('time step')

% figure(4)
% plot(c./fmax*1e9)
% ylabel('Modulated wavelength (nm)')
% xlabel('time step')

close all
Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector
S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);
X = S + 2*randn(size(t));
plot(1000*t(1:50),X(1:50))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (milliseconds)')
ylabel('X(t)')
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
