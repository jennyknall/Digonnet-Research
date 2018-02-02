function [N2,N1,N2j,N1i] = popInversionVsWavelength(P_l,P_c,cs_abs,cs_ems)

global E11 E12 E13 E14 E21 E22 E23 E1 E2 g2 g1 tau21 N0 Area len;

%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%
h = 6.63e-34; %J*s
c = 3e8;
kT = 4.11e-21; %J

%frequency of the laser/cooling pumps
freq_l = c./cs_abs(:,1);
freq_c = [c/990e-9]; %c./cs_abs(:,1); %CHANGE: sweep only over a couple cooling pump wavelengths

%%%%% Calculating N2 and N1 %%%%%%%%%%%
L = P_l/Area/h./freq_l;
C = P_c/Area/h./freq_c;

N2 = cell(1,length(freq_c));
N1 = cell(1,length(freq_c));
for i = 1:length(freq_c) %for each cooling pump frequency
    CFI = 91; %%%find index of cs_abs that cooresponds to freq_c(i)
    N2{i} = (L.*cs_abs(:,2)+C(i)*cs_abs(CFI,2))./(L.*(cs_abs(:,2)+cs_ems(:,2))+C(i)*(cs_abs(CFI,2)+cs_ems(CFI,2))+1/tau21)*N0;
    N1{i} = N0-N2{i};
end

%%%%%% Calculating N2j and N1i %%%%%%%%%
N2j = cell(g2,length(freq_c)); %dim-> height:population of each sublevel for each laser pumping power
                            %      length:different cooling pumping power for each collumn
N1i = cell(g1,length(freq_c));

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

for p = 1:length(freq_c) %for each cooling pump power, find the populations of each sublevel.
    for j = 1:g2
        delE = (E2(j)-E21)*h*c;
        N2j{j,p} = exp(-delE/kT)*N2{p}/N2_N21;
    end
    for i = 1:g1
        delE = (E1(i)-E11)*h*c;
        N1i{i,p} = exp(-delE/kT)*N1{p}/N1_N11;
    end
end

%%%%%%%% GRAPHS %%%%%%%%%%%%%
%N1
figure(1)
hold on
grid on
for i = 1:length(N1)
    plot(c./freq_l*1e9,N1{i}/N0,'DisplayName',sprintf('lambda_c = %g nm',c/freq_c(i)*1e9))
end
xlabel('Laser pump wavelength (nm)');
ylabel('N1/N0');
title(sprintf('N1 vs. Laser pumping wavelength for P_l  = %g W, P_c = %g W',P_l,P_c))
legend('show')

%N2
figure(2)
hold on
grid on
for i = 1:length(N2)
    plot(c./freq_l*1e9,N2{i}/N0,'DisplayName',sprintf('lambda_c = %g nm',c/freq_c(i)*1e9))
end
xlabel('Laser pump wavelength (nm)');
ylabel('N2/N0');
title(sprintf('N2 vs. Laser pumping wavelength for P_l  = %g W, P_c = %g W',P_l,P_c))
legend('show')

%N2-N1
figure(3)
hold on
grid on
for i = 1:length(N2)
    plot(c./freq_l*1e9,(N2{i}-N1{i})/N0,'DisplayName',sprintf('lambda_c = %g nm',c/freq_c(i)*1e9))
end
xlabel('Laser pump wavelength (nm)');
ylabel('(N2-N1)/N0');
title(sprintf('N2-N1 vs. Laser pumping wavelength for P_l  = %g W, P_c = %g W',P_l,P_c))
legend('show','Location','southeast')

%N21-N11
figure(4)
hold on
grid on
for i = 1:length(freq_c) %for each cooling pump power
    plot(c./freq_l*1e9,(N2j{1,i}-N1i{1,i})/N0,'DisplayName',sprintf('lambda_c = %g nm',c/freq_c(i)*1e9))
end
xlabel('Laser pump wavelength (nm)');
ylabel('(N21-N11)/N0');
title(sprintf('N21-N11 vs. Laser pumping wavelength for P_l  = %g W, P_c = %g W',P_l,P_c))
legend('show','Location','southeast')

%N2j-N1i
figure(5)
hold on
grid on
%ylim([-0.2,0.6]);
CPI= 1; %Cooling Power Index
numGraphs = 1;
for i = 1:g1
    for j = 1:g2
        plot(c./freq_l*1e9,(N2j{j,CPI}-N1i{i,CPI})/N0,'lineWidth',ceil(numGraphs/6)*2,'DisplayName',sprintf('N2%d - N1%d',j,i));
        numGraphs = numGraphs+1;
    end
end
xlabel('Laser Pump wavelength (nm)');
ylabel('(N2j-N1i)/N0');
title(sprintf('N2j-N1i vs. Laser pumping wavelength, lamda_c = %g nm for P_l  = %g W, P_c = %g W',c/freq_c(CPI)*1e9,P_l,P_c))
legend('show','Location','southeast');


end