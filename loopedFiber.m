%calculate heat extraction for looped fiber case
%r_looped = determineRadius(loops,b);
%h_looped = 40; %W/m^2/K
zNew = linspace(0,len,length(z)*100);
dQz_dtInterp = interp1(z,dQz_dt2,zNew);
len_new = .92;
index = floor(len_new/len*length(zNew));
dQ_new = dQz_dtInterp(1,1:index);
loops = 4;
segment = floor(length(dQ_new)/loops);
dQ_looped = zeros(1,segment);

for i = 1:loops
    seg = (i-1)*segment+1:i*segment;
    dQ_looped = dQ_looped+dQ_new(seg);
end
lastBit = loops*segment+1:length(dQ_new);
zeroarray = zeros(1,segment-length(lastBit));
dQ_looped = dQ_looped+[dQ_new(lastBit) zeroarray];
dQcheck = trapz(zNew(1:segment),dQ_looped)
%newTempSeg = dQ_looped/2/pi/r_looped/h_looped;

%dQ as a function of z for loop
figure(7)
grid on 
hold on
plot(zNew(1:segment),dQ_looped,'DisplayName',sprintf('%g mW',Pp0(p)*1e3));
plot(zNew(1:index),dQ_new,'DisplayName',sprintf('%g mW',Pp0(p)*1e3));
xlabel('z (m)');
ylabel('Heat extracted per unit length (W/m)');
title('Heat extraction along the fiber for looped case');
