concentration = 0.2; 

densityIon = 2.39; %g/cm^3
densityHost = 2.65; %g/cm^3

numAtoms = 2;

massIon = 2*30.974+5*15.999; %g/mol

V1 = 1; %cm^3
V2 = (densityIon*V1-concentration*densityIon*V1)/densityHost/concentration;

Vtot = V1+V2;

AVN = 6.022e23; %avagadros number

N0 = densityIon*V1/Vtot * 1/massIon * AVN * numAtoms * 1e6 %m^-3

