function [N2,N1,N2j,N1i] = popInversionVsPower(lambda_l,lambda_c,cs_al,cs_el,cs_ac,cs_ec,P_l,P_c)

global E11 E12 E13 E14 E21 E22 E23 E1 E2 g2 g1 tau21 N0 Area len;

%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%
h = 6.63e-34; %J*s
c = 3e8;
kT = 4.11e-21; %J

freq_l = c/lambda_l;
freq_c = c/lambda_c;

%%%%% Calculating N2 and N1 %%%%%%%%%%%
L = P_l/Area/h/freq_l;
C = P_c/Area/h/freq_c;

N2 = cell(1,length(P_c));
N1 = cell(1,length(P_c));
for i = 1:length(P_c)
    N2{i} = (L*cs_al+C(i)*cs_ac)./(L*(cs_al+cs_el)+C(i)*(cs_ac+cs_ec)+1/tau21)*N0;
    N1{i} = N0-N2{i};
end

%%%%%% Calculating N2j and N1i %%%%%%%%%
N2j = cell(g2,length(P_c)); %dim-> height:population of each sublevel for each laser pumping power
                            %      length:different cooling pumping power for each collumn
N1i = cell(g1,length(P_c));

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

for p = 1:length(P_c) %for each cooling pump power, find the populations of each sublevel.
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
    plot(P_l,N1{i}/N0,'DisplayName',sprintf('P_c = %f W',P_c(i)))
end
xlabel('Laser pump power (W)');
ylabel('N1/N0');
title('N1 vs. Laser pumping power')
legend('show')

%N2
figure(2)
hold on
grid on
for i = 1:length(N2)
    plot(P_l,N2{i}/N0,'DisplayName',sprintf('P_c = %f W',P_c(i)))
end
xlabel('Laser pump power (W)');
ylabel('N2/N0');
title('N2 vs. Laser pumping power')
legend('show','Location','southeast')

%N2-N1
figure(3)
hold on
grid on
for i = 1:length(N2)
    plot(P_l,(N2{i}-N1{i})/N0,'DisplayName',sprintf('P_c = %f W',P_c(i)))
end
xlabel('Laser pump power (W)');
ylabel('(N2-N1)/N0');
title('N2-N1 vs. Laser pumping power')
legend('show','Location','southeast')

%N21-N11
figure(4)
hold on
grid on
for i = 1:length(P_c) %for each cooling pump power
    plot(P_l,(N2j{1,i}-N1i{1,i})/N0,'DisplayName',sprintf('P_c = %g W',P_c(i)))
end
xlabel('Laser pump power (W)');
ylabel('(N21-N11)/N0');
title('N21-N11 vs. Laser pumping power')
legend('show','Location','southeast')

%N2j-N1i
figure(5)
hold on
grid on
ylim([-0.2,0.6]);
CPI= 5; %Cooling Power Index
numGraphs = 1;
for i = 1:g1
    for j = 1:g2
        plot(P_l,(N2j{j,CPI}-N1i{i,CPI})/N0,'lineWidth',ceil(numGraphs/6)*2,'DisplayName',sprintf('N2%d - N1%d',j,i));
        numGraphs = numGraphs+1;
    end
end
xlabel('Laser Pump Power (W)');
ylabel('(N2j-N1i)/N0');
title(sprintf('N2j-N1i vs. Laser pumping power, P_c = %g W',P_c(CPI)))
legend('show','Location','southeast');

%N2j-N1i
figure(6)
hold on
grid on
ylim([-0.2,0.6]);
CPI= 1; %Cooling Power Index
numGraphs = 1;
for i = 1:g1
    for j = 1:g2
        plot(P_l,(N2j{j,CPI}-N1i{i,CPI})/N0,'lineWidth',ceil(numGraphs/6)*2,'DisplayName',sprintf('N2%d - N1%d',j,i));
        numGraphs = numGraphs+1;
    end
end
xlabel('Laser Pump Power (W)');
ylabel('(N2j-N1i)/N0');
title(sprintf('N2j-N1i vs. Laser pumping power, P_c = %g W',P_c(CPI)))
legend('show','Location','southeast');





end