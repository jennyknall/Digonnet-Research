function [N2,N1,N2j,N1i] = popInversionVsWavelength3(Pl0,Pc0,lambda_f,cs_abs,cs_ems)
%%This function takes into account the power change along the fiber for TWO PUMPS
%%Also calculates Q vs. Wavelength

global E11 E12 E13 E14 E21 E22 E23 E1 E2 g2 g1 tau21 N0 Area len core_radius loss;

%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%
h = 6.63e-34; %J*s
c = 3e8;
kT = 4.11e-21; %J
f = 1/Area;

%frequency of the lasing and cooling pumps
freq_l = c./cs_abs(:,1);
lambda_c = 1005e-9;
freq_c = c/lambda_c;
cIndex = 136; %index of cs_abs/cs_ems that cooresponds to lambda_c

%%%%% Calculating N2 and N1 %%%%%%%%%%%
N2 = NaN(1,length(freq_l));
N1 = NaN(1,length(freq_l));
Qtotal = NaN(1,length(freq_l));
Q1 = NaN(1,length(freq_l));
Q2 = NaN(1,length(freq_l));
dz = 0.5; %change in power calculated in increments of dz
z = 0:dz:len;
for i = 1:length(freq_l)
    
    %%caluclate P(z)
    [Plz, Pcz] = changeInPumpPower2(Pl0,Pc0,cs_abs(i,1),lambda_c,cs_abs(i,2),cs_ems(i,2),cs_abs(cIndex,2),cs_ems(cIndex,2),dz); 
    
    Isatl = h*freq_l(i)/(cs_abs(i,2)+cs_ems(i,2))/tau21;
    Isatc = h*freq_c/(cs_abs(cIndex,2)+cs_ems(cIndex,2))/tau21;
    
    %calculate N2(z), N1(z)
    N2z = (Plz*f/h/freq_l(i)*cs_abs(i,2)+Pcz*f/h/freq_c*cs_abs(cIndex,2))*tau21./(1+Plz*f/Isatl+Pcz*f/Isatc)*N0;
    N1z = N0-N2z;
    
    %calculate Qz
    dQ_dV = Plz*f.*(cs_abs(i,2)*N1z-cs_ems(i,2)*N2z+loss)+Pcz*f.*(cs_abs(cIndex,2)*N1z-cs_ems(cIndex,2)*N2z+loss)-N2z/tau21*h*c/lambda_f;
    dQ1_dV = Plz*f.*(cs_abs(i,2)*N1z-cs_ems(i,2)*N2z+loss)+Pcz*f.*(cs_abs(cIndex,2)*N1z-cs_ems(cIndex,2)*N2z+loss);
    dQ2_dV = -N2z/tau21*h*c/lambda_f;
    Q1(i) = pi*core_radius^2*trapz(z,dQ1_dV);
    Q2(i) = pi*core_radius^2*trapz(z,dQ2_dV);
    Qtotal(i) = pi*core_radius^2*trapz(z,dQ_dV);
    
    N2(i) = trapz(z,N2z)/N0/len;
    N1(i) = trapz(z,N1z)/N0/len;
end


%%%%%% Calculating N2j and N1i %%%%%%%%%
N2j = zeros(g2,length(N2)); %dim-> height:population of each sublevel for each wavelength
                            %      length:different wavelengths
N1i = zeros(g1,length(N1));

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

%find population for each sublevel
for j = 1:g2
    delE = (E2(j)-E21)*h*c;
    N2j(j,:) = exp(-delE/kT)*N2/N2_N21;
end
for i = 1:g1
    delE = (E1(i)-E11)*h*c;
    N1i(i,:) = exp(-delE/kT)*N1/N1_N11;
end

%normalizing Q to 1 so that it can be graphed on the same graphs as N1, N2
Qnorm = Qtotal/max(Qtotal);

%%%%%%% GRAPHS %%%%%%%%%%%%%
%N1
figure(1)
grid on
plot(c./freq_l*1e9,N1)
xlabel('Laser pump wavelength (nm)');
ylabel('N1/N0');
title(sprintf('N1 vs. Laser pumping wavelength for Pl0  = %g W, Pc0 = %g',Pl0,Pc0))

%N2
figure(2)
grid on
plot(c./freq_l*1e9,N2)
xlabel('Laser pump wavelength (nm)');
ylabel('N2/N0');
title(sprintf('N2 vs. Laser pumping wavelength for Pl0  = %g W, Pc0 = %g',Pl0,Pc0))

%N2-N1
figure(3)
grid on
hold on
plot(c./freq_l*1e9,(N2-N1))
plot(c./freq_l*1e9,Qtotal)
xlabel('Laser pump wavelength (nm)');
ylabel('(N2-N1)/N0 and Q');
title(sprintf('N2-N1 vs. Laser pumping wavelength for Pl0  = %g W, Pc0 = %g',Pl0,Pc0))
legend('(N2-N1)/N0','Q')

%N21-N11
figure(4)
grid on
hold on
plot(c./freq_l*1e9,(N2j(1,:)-N1i(1,:)))
plot(c./freq_l*1e9,Qtotal)
xlabel('Laser pump wavelength (nm)');
ylabel('(N21-N11)/N0 and Q');
title(sprintf('N21-N11 vs. Laser pumping wavelength for Pl0  = %g W, Pc0 = %g',Pl0,Pc0))
legend('(N21-N11)/N0','Q')

%N2j-N1i
figure(5)
hold on
grid on
%ylim([-0.2,0.6]);
numGraphs = 1;
for i = 1:g1
    for j = 1:g2
        plot(c./freq_l*1e9,(N2j(j,:)-N1i(i,:)),'lineWidth',ceil(numGraphs/6)*2,'DisplayName',sprintf('N2%d - N1%d',j,i));
        numGraphs = numGraphs+1;
    end
end
plot(c./freq_l*1e9,Qtotal,'k','lineWidth',4,'DisplayName','Q')
xlabel('Laser Pump wavelength (nm)');
ylabel('(N2j-N1i)/N0 and Q');
title(sprintf('N2j-N1i vs. Laser pumping wavelength for Pl0  = %g W, Pc0 = %g',Pl0,Pc0))
legend('show','Location','southeast');

figure(6)
grid on
plot(c./freq_l*1e9,Qtotal)
xlabel('Laser pump wavelength (nm)');
ylabel('Q');
title(sprintf('Q vs. Laser pumping wavelength for Pl0  = %g W, Pc0 = %g',Pl0,Pc0))

figure(7)
grid on
plot(c./freq_l*1e9,Q1)
xlabel('Laser pump wavelength (nm)');
ylabel('Q1');
title(sprintf('Q1 vs. Laser pumping wavelength for Pl0  = %g W, Pc0 = %g',Pl0,Pc0))

figure(8)
grid on
plot(c./freq_l*1e9,Q2)
xlabel('Laser pump wavelength (nm)');
ylabel('Q2');
title(sprintf('Q2 vs. Laser pumping wavelength for Pl0  = %g W, Pc0 = %g',Pl0,Pc0))


end