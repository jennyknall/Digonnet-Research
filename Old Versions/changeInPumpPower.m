function P = changeInPumpPower(P0,lambda,cs_a,cs_e,dz)
%%This function takes into account the power change along the fiber for ONE PUMP


global tau21 N0 Area len loss;

%Constants
h = 6.63e-34; %J*s
c = 3e8;

%frequency of laser/cooling pump
freq = c/lambda;

Isat = h*freq/(cs_a+cs_e)/tau21;
f = 1/Area;
gamma = 1;
                                     
z = 0:dz:len;                                         
P = zeros(1,length(z)); 
P(1) = P0; %W                                      % initial pump power
dP_dz = @(z,p) ((cs_a+cs_e)*(p*f/h/freq*cs_a*tau21)/(1+p*f/Isat)*N0-cs_a*N0)*p*gamma-loss*p*gamma;  

for i=1:(length(z)-1)                              % calculation loop
    k_1 = dP_dz(z(i),P(i));
    k_2 = dP_dz(z(i)+0.5*dz,P(i)+0.5*dz*k_1);
    k_3 = dP_dz((z(i)+0.5*dz),(P(i)+0.5*dz*k_2));
    k_4 = dP_dz((z(i)+dz),(P(i)+k_3*dz));

    P(i+1) = P(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dz;  % main equation
end

%graphs
plot(z,P);
grid on
xlabel('z (m)');
ylabel('P(z) (W)');
title('Change in Power along fiber');