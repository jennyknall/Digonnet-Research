load('h_values.mat')

b = 62.5e-6;
thermalC = 81.4;

delT = 29.9;

delQ = delT*2*pi*b*thermalC;

thermalC = hINTERP(find(abs(hINTERP(:,1) - delT)-1e-9 < 0),2);
delT = delQ/2/pi/b/thermalC

thermalC = hINTERP(find(abs(hINTERP(:,1) - round(delT,1))-1e-9 < 0),2);
delT = delQ/2/pi/b/thermalC

thermalC = hINTERP(find(abs(hINTERP(:,1) - round(delT,1))-1e-9 < 0),2);
delT = delQ/2/pi/b/thermalC