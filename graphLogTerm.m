loss = maxLoss/1000:maxLoss/1000:maxLoss*etta;
lossdB = loss*10/log(10);
pump = optPump(loss/etta);

exactLog = log((1+pump/PsatP*(1-etta))./(1+pump/PsatP));
approxLog = log(1-etta)*ones(1,length(exactLog));

figure(1)
hold on
plot(pump,exactLog)
plot(pump,approxLog)
xlabel('Pin = Popt (W)')
legend('exact log term','log(1-etta)')

figure(2)
hold on
plot(lossdB,exactLog)
plot(lossdB,approxLog)
xlabel('loss (dB/m)');
legend('exact log term','log(1-etta)')

figure(3)
plot(lossdB,pump/PsatP)
xlabel('loss (dB/m)')
ylabel('Pin/PsatP')
figure(2)

ratio = 1:100;
log1 = log((1+ratio*(1-etta))./(1+ratio));
log2 = log(1-etta)*ones(1,length(log1));
figure(4)
hold on
plot(ratio,log1);
plot(ratio,log2);
xlabel('Pin/Psat');
legend('exact log term','log(1-etta)')

ratio = 1:100;
log1 = (1+ratio*(1-etta))./(1+ratio);
log2 = (1-etta)*ones(1,length(log1));
figure(5)
hold on
plot(ratio,log1);
plot(ratio,log2);
xlabel('Pin/Psat');
legend('exact log term','1-etta')

