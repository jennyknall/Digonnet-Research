function [N2,N1,N2j,N1i] = popInversionVsPower3(lambda_l,lambda_c,lambda_f,cs_al,cs_el,cs_ac,cs_ec,Pl0,Pc0)
%%This function takes into account the power change along the fiber for TWO PUMPS
%%Also calculates Q vs. Wavelength

global E11 E12 E13 E14 E21 E22 E23 E1 E2 g2 g1 tau21 N0 Area len core_radius loss;

%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%
h = 6.63e-34; %J*s
c = 3e8;
kT = 4.11e-21; %J
f = 1/Area;

freq_l = c/lambda_l;
freq_c = c/lambda_c;

%%%%% Calculating N2 and N1 %%%%%%%%%%%
N2 = cell(1,length(Pc0));
N1 = cell(1,length(Pc0));
Q = cell(1,length(Pc0));
dz = 0.5;  %change in power calculated in increments of dz
z = 0:dz:len;
for i = 1:length(Pc0) %for each initial cooling pump power
    N2{i} = NaN(1,length(Pl0));
    Q{i} = NaN(1,length(Pl0));
    for j = 1:length(Pl0) %for each initial lasing pump power
        %calculate Pl(z),Pc(z)
        [Plz, Pcz] = changeInPumpPower2(Pl0(j),Pc0(i),lambda_l,lambda_c,cs_al,cs_el,cs_ac,cs_ec,dz); 
        
        Isatl = h*freq_l/(cs_al+cs_el)/tau21;
        Isatc = h*freq_c/(cs_ac+cs_ec)/tau21;
        
        %calculate N2(z), N1(z)
        N2z = (Plz*f/h/freq_l*cs_al+Pcz*f/h/freq_c*cs_ac)*tau21./(1+Plz*f/Isatl+Pcz*f/Isatc)*N0;
        N1z = N0-N2z;
        
        %calculate Qz
        dQ_dV = Plz*f.*(cs_al*N1z-cs_el*N2z+loss)+Pcz*f.*(cs_ac*N1z-cs_ec*N2z+loss)-N2z/tau21*h*c/lambda_f;
        Q{i}(j) = pi*core_radius^2*trapz(z,dQ_dV);
        
        N2{i}(j) = trapz(z,N2z)/N0/len;
        N1{i}(j) = trapz(z,N1z)/N0/len;
    end
end

%%%%%% Calculating N2j and N1i %%%%%%%%%
N2j = cell(g2,length(Pc0)); %dim-> height:population of each sublevel for each initial laser pumping power
                            %      length:different initial cooling pumping power for each collumn
N1i = cell(g1,length(Pc0));

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

for p = 1:length(Pc0) %for each cooling pump power, find the populations of each sublevel.
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
    plot(Pl0,N1{i},'DisplayName',sprintf('Pc0 = %g W',Pc0(i)))
    plot(Pl0,Q{i},'DisplayName',sprintf('Q for Pc0 = %g W',Pc0(i)),'lineWidth',4)
end
xlabel('Initial laser pump power (W)');
ylabel('N1/N0 and Q');
title('N1 vs. Laser pumping power')
legend('show')

%N2
figure(2)
hold on
grid on
for i = 1:length(N2)
    plot(Pl0,N2{i},'DisplayName',sprintf('Pc0 = %g W',Pc0(i)))
    plot(Pl0,Q{i},'DisplayName',sprintf('Q for Pc0 = %g W',Pc0(i)),'lineWidth',4)
end
xlabel('Initial laser pump power (W)');
ylabel('N2/N0 and Q');
title('N2 vs. Laser pumping power')
legend('show','Location','southeast')

%N2-N1
figure(3)
hold on
grid on
for i = 1:length(N2)
    plot(Pl0,N2{i}-N1{i},'DisplayName',sprintf('Pc0 = %g W',Pc0(i)))
    plot(Pl0,Q{i},'DisplayName',sprintf('Q for Pc0 = %g W',Pc0(i)),'lineWidth',4)
end
xlabel('Initial laser pump power (W)');
ylabel('(N2-N1)/N0 and Q');
title('N2-N1 vs. Laser pumping power')
legend('show','Location','southeast')

%N21-N11
figure(4)
hold on
grid on
for i = 1:length(Pc0) %for each cooling pump power
    plot(Pl0,N2j{1,i}-N1i{1,i},'DisplayName',sprintf('Pc0 = %g W',Pc0(i)))
    plot(Pl0,Q{i},'DisplayName',sprintf('Q for Pc0 = %g W',Pc0(i)),'lineWidth',4)
end
xlabel('Initial laser pump power (W)');
ylabel('(N21-N11)/N0 and Q');
title('N21-N11 vs. Laser pumping power')
legend('show','Location','southeast')

%N2j-N1i
figure(5)
hold on
grid on
ylim([-0.2,0.6]);
CPI= 1; %Cooling Power Index
numGraphs = 1;
for i = 1:g1
    for j = 1:g2
        plot(Pl0,N2j{j,CPI}-N1i{i,CPI},'lineWidth',ceil(numGraphs/6)*2,'DisplayName',sprintf('N2%d - N1%d',j,i));
        numGraphs = numGraphs+1;
    end
end
plot(Pl0,Q{CPI},'k','DisplayName',sprintf('Q for Pc0 = %g W',Pc0(CPI)),'lineWidth',4)
xlabel('Initial laser Pump Power (W)');
ylabel('(N2j-N1i)/N0');
title(sprintf('N2j-N1i vs. Laser pumping power, Pc0 = %g W',Pc0(CPI)))
legend('show','Location','southeast');

figure(6)
hold on
grid on
for i = 1:length(N1)
    plot(Pl0,Q{i},'DisplayName',sprintf('Q for Pc0 = %g W',Pc0(i)))
end
xlabel('Initial laser pump power (W)');
ylabel('Q');
title(sprintf('Q vs. Laser pumping power, pump wavelength = %g nm',lambda_l*1e9))
legend('show')

end