% need to run plotLoss90vsN0.m before this script inorder for all variables to be defined

close all

approxFac = .948;

D = percentCooling*2*core_radius^2/w^2+approxFac*log(1-etta);
%Kplus = 1/2+(2-etta)/2/etta*D+1/2*sqrt(1+2*(2-etta)/etta*D+D.^2);
K = 1/2+(2-etta)/2/etta*D-1/2*sqrt(1+2*(2-etta)/etta*D+D.^2);

changeD = (D(1)-D(end))/max(abs(D))
changeK = (K(1)-K(end))/max(abs(K))

figure(1)
plot(approxFac,D)
ylabel('D')

figure(2)
plot(approxFac,K)
ylabel('Kminus');
