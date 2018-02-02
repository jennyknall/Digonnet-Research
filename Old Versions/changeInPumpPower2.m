function [Pl,Pc] = changeInPumpPower2(Pl0,Pc0,lambda_l,lambda_c,cs_al,cs_el,cs_ac,cs_ec,dz)
%%This function takes into account the power change along the fiber for TWO PUMPS

global tau21 N0 Area len loss;

%Constants
h = 6.63e-34; %J*s
c = 3e8;

%frequency of laser/cooling pump
freq_l = c/lambda_l;
freq_c = c/lambda_c;

Isatl = h*freq_l/(cs_al+cs_el)/tau21;
Isatc = h*freq_c/(cs_ac+cs_ec)/tau21;
f = 1/Area;
gamma = 1;
                                     
z = 0:dz:len;                                         
Pl = zeros(1,length(z)); 
Pc = zeros(1,length(z));
Pl(1) = Pl0; %W                                      % initial laser pump power
Pc(1) = Pc0; %W                                      % initial cooling pump power
%laser pumping power
dPl_dz = @(z,pl,pc) ((cs_al+cs_el)*(pl*f/h/freq_l*cs_al+pc*f/h/freq_c*cs_ac)...
                        *tau21/(1+pl*f/Isatl+pc*f/Isatc)*N0-cs_al*N0)*pl*gamma-loss*pl*gamma;  
%cooling pumping power
dPc_dz = @(z,pl,pc) ((cs_ac+cs_ec)*(pl*f/h/freq_l*cs_al+pc*f/h/freq_c*cs_ac)...
                        *tau21/(1+pl*f/Isatl+pc*f/Isatc)*N0-cs_ac*N0)*pc*gamma-loss*pc*gamma; 
for i=1:(length(z)-1)                              % calculation loop
    k_1 = dPl_dz(z(i),Pl(i),Pc(i));
    j_1 = dPc_dz(z(i),Pl(i),Pc(i));
    k_2 = dPl_dz(z(i)+0.5*dz,Pl(i)+0.5*dz*k_1,Pc(i)+0.5*dz*j_1);
    j_2 = dPc_dz(z(i)+0.5*dz,Pl(i)+0.5*dz*k_1,Pc(i)+0.5*dz*j_1);
    k_3 = dPl_dz((z(i)+0.5*dz),(Pl(i)+0.5*dz*k_2),(Pc(i)+0.5*dz*j_2));
    j_3 = dPc_dz((z(i)+0.5*dz),(Pl(i)+0.5*dz*k_2),(Pc(i)+0.5*dz*j_2));
    k_4 = dPl_dz((z(i)+dz),(Pl(i)+k_3*dz),(Pc(i)+j_3*dz));
    j_4 = dPc_dz((z(i)+dz),(Pl(i)+k_3*dz),(Pc(i)+j_3*dz));

    Pl(i+1) = Pl(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dz;  % main equation
    Pc(i+1) = Pc(i) + (1/6)*(j_1+2*j_2+2*j_3+j_4)*dz;  % main equation
end

%graphs
% 
% figure(7)
% plot(z,Pc);
% xlabel('z (m)');
% ylabel('Pc(z) (W)');
% title('Change in cooling pump Power along fiber');

% figure(6)
% plot(z,Pl);
% xlabel('z (m)');
% ylabel('Pl(z) (W)');
% title('Change in laser pump Power along fiber');
