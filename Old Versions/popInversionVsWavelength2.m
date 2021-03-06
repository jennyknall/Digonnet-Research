function [N2,N1,N2j,N1i,Q] = popInversionVsWavelength2(P0,lambda_f,cs_abs,cs_ems)
%%This function takes into account the power change along the fiber for ONE PUMP
%%Also calculates Q vs. Wavelength. APPROXIMATES intensity in
%fiber as P(z)/Area

global E11 E21 E1 E2 g2 g1 tau21 N0 Area len core_radius loss;

%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%
h = 6.63e-34; %J*s
c = 3e8;
kT = 4.11e-21; %J

%frequency of the pump
freq = c./cs_abs(:,1);


%%%%% Calculating N2 and N1 %%%%%%%%%%%
N2 = NaN(1,length(freq));
N1 = NaN(1,length(freq));
Q = NaN(1,length(freq));
for i = 1:length(freq)
    
    %%caluclate P(z)
    f = 1/Area;
    dz = 0.5; %change in power calculated in increments of dz
    z = 0:dz:len;
    Pz = changeInPumpPower(P0,cs_abs(i,1),cs_abs(i,2),cs_ems(i,2),dz); 
    
    Isat = h*freq(i)/(cs_abs(i,2)+cs_ems(i,2))/tau21;
    
    %calculate N2(z), N1(z)
    N2z = Pz*f/h/freq(i)*cs_abs(i,2)*tau21./(1+Pz*f/Isat)*N0;
    N1z = N0-N2z;
    
    %calculate Qz
    dQ_dV = Pz*f.*(cs_abs(i,2)*N1z-cs_ems(i,2)*N2z+loss)-N2z/tau21*h*c/lambda_f;
    Q(i) = pi*core_radius^2*trapz(z,dQ_dV);

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

for j = 1:g2
    delE = (E2(j)-E21)*h*c;
    N2j(j,:) = exp(-delE/kT)*N2/N2_N21;
end
for i = 1:g1
    delE = (E1(i)-E11)*h*c;
    N1i(i,:) = exp(-delE/kT)*N1/N1_N11;
end

%normalizing Q to 1 so that it can be graphed on the same graphs as N1, N2
Qnorm = Q/max(Q);

%%%%%%% GRAPHS %%%%%%%%%%%%%
%N1
figure(1)
grid on
plot(c./freq*1e9,N1)
xlabel('Laser pump wavelength (nm)');
ylabel('N1/N0');
title(sprintf('N1 vs. Laser pumping wavelength for P0  = %g W',P0))

%N2
figure(2)
grid on
hold on
plot(c./freq*1e9,N2)
plot(c./freq*1e9,N1)
xlabel('Laser pump wavelength (nm)');
ylabel('N2/N0');
title(sprintf('N2 vs. Laser pumping wavelength for P0  = %g W',P0))

%N2-N1 and Q
figure(3)
grid on
hold on
plot(c./freq*1e9,(N2-N1))
plot(c./freq*1e9,Qnorm)
xlabel('Laser pump wavelength (nm)');
ylabel('(N2-N1)/N0 and Q');
title(sprintf('N2-N1 vs. Laser pumping wavelength for P0  = %g W',P0))
legend('N2-N1','normalized Q')

%N21-N11
figure(4)
grid on
plot(c./freq*1e9,(N2j(1,:)-N1i(1,:)))
xlabel('Laser pump wavelength (nm)');
ylabel('(N21-N11)/N0');
title(sprintf('N21-N11 vs. Laser pumping wavelength for P0  = %g W',P0))

%N2j-N1i and Q
figure(5)
hold on
grid on
%ylim([-0.2,0.6]);
numGraphs = 1;
for i = 1:g1
    for j = 1:g2
        plot(c./freq*1e9,(N2j(j,:)-N1i(i,:)),'lineWidth',ceil(numGraphs/6)*2,'DisplayName',sprintf('N2%d - N1%d',j,i));
        numGraphs = numGraphs+1;
    end
end
plot(c./freq*1e9,Qnorm,'k','lineWidth',4,'DisplayName','Qnorm')
xlabel('Laser Pump wavelength (nm)');
ylabel('(N2j-N1i)/N0');
title(sprintf('N2j-N1i vs. Laser pumping wavelength for P0  = %g W',P0))
legend('show','Location','southeast');

%Q
figure(6)
grid on
plot(c./freq*1e9,Q)
xlabel('Laser pump wavelength (nm)');
ylabel('Q (W)');
title(sprintf('Q vs. Laser pumping wavelength for P0  = %g W',P0))


end