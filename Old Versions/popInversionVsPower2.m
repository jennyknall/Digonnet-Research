function [N2,N1,N2j,N1i] = popInversionVsPower2(lambda,cs_abs,cs_ems,P0)
%%This function takes into account the power change along the fiber for ONE PUMPS

global E11 E21 E1 E2 g2 g1 tau21 N0 Area len;

%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%
h = 6.63e-34; %J*s
c = 3e8;
kT = 4.11e-21; %J

freq = c/lambda;

%%%%% Calculating N2 and N1 %%%%%%%%%%%
N2 = NaN(1,length(P0));
N1 = NaN(1,length(P0));
dz = 0.5;  %change in power calculated in increments of dz
z = 0:dz:len;
for j = 1:length(P0) %for each initial lasing pump power
    %calculate Pl(z)
    Pz = changeInPumpPower(P0(j),lambda,cs_abs,cs_ems,dz); 

    f = 1/Area;
    Isat = h*freq/(cs_abs+cs_ems)/tau21;

    %calculate N2(z), N1(z)
    N2z = Pz*f/h/freq*cs_abs*tau21./(1+Pz*f/Isat)*N0;
    N1z = N0-N2z;

    N2(j) = trapz(z,N2z)/N0/len;
    N1(j) = trapz(z,N1z)/N0/len;
end

%%%%%% Calculating N2j and N1i %%%%%%%%%
N2j = NaN(g2,length(N2)); %dim-> height:population of each sublevel for each power
                            %      length:different powers
N1i = NaN(g1,length(N1));

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


%%%%%%%% GRAPHS %%%%%%%%%%%%%
%N1
figure(1)
grid on
plot(P0,N1)
xlabel('Initial laser pump power (W)');
ylabel('N1/N0');
title('N1 vs. Laser pumping power')

%N2
figure(2)
grid on
plot(P0,N2)
xlabel('Initial laser pump power (W)');
ylabel('N2/N0');
title('N2 vs. Laser pumping power')

%N2-N1
figure(3)
grid on
plot(P0,N2-N1)
xlabel('Initial laser pump power (W)');
ylabel('(N2-N1)/N0');
title('N2-N1 vs. Laser pumping power')

%N21-N11
figure(4)
grid on
plot(P0,N2j(1,:)-N1i(1,:))
xlabel('Initial laser pump power (W)');
ylabel('(N21-N11)/N0');
title('N21-N11 vs. Laser pumping power')

%N2j-N1i
figure(5)
hold on
grid on
ylim([-0.2,0.6]);
numGraphs = 1;
for i = 1:g1
    for j = 1:g2
        plot(P0,N2j(j,:)-N1i(i,:),'lineWidth',ceil(numGraphs/6)*2,'DisplayName',sprintf('N2%d - N1%d',j,i));
        numGraphs = numGraphs+1;
    end
end
xlabel('Initial laser Pump Power (W)');
ylabel('(N2j-N1i)/N0');
title('N2j-N1i vs. Laser pumping power')
legend('show','Location','southeast');



end