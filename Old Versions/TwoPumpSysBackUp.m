clear all;
close all;

%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%
h = 6.63e-34; %J*s
c = 3e8;
kT = 4.11e-21; %J

%energy levels, m^-1
E11 = 0;
E12 = 338e2;
E13 = 445e2;
E14 = 872e2;
E21 = 10288e2;
E22 = 10438e2;
E23 = 11038e2;

E1 = [E11 E12 E13 E14];
E2 = [E21 E22 E23];

g2 = 3;
g1 = 4;

%absorption/emission cross sectional areas
cs_abs = xlsread('absorption_00.xlsx');
cs_ems = xlsread('emission_00.xlsx');
cs_al = 23e-26; %m^2 
cs_el = 1e-26; %m^2 
cs_ac = 26e-26; %m^2 
cs_ec = 6e-26; %m^2 
tau21 = 1e-3; %s
N0 = 1.1e25; %m^-3;
Area = 1e-10; %m^2, pump area
L = 10; %m; length of fiber

%%%%%%%%%%%% Sweeping Variables %%%%%%%%%%%%%%%%%%%
%laser/cooling pump power
P_l = 0.001:0.001:2; %W
P_c = 0.1:0.1:0.5; %W

%waveleghth of laser/cooling pump
lambda_l = 1/(E23-E12); %1/(E23-E11):e-9:1/(E22-E11);
freq_l = c/lambda_l;
lambda_c = 1/(E21-E12); %1/(E22-E11):e-9:1/(E21-E11);
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
for i = 1:length(N2j) %for each cooling pump power
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

