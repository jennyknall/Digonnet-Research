clear all; close all;

syms r;
NA = 0.2;
w = 2.2e-6;
lam = 1060e-9;
r = (1:.01:2)*1e-6;
y = r.*(0.65+1.619./(2*pi*r/lam*NA).^1.5+2.879./(2*pi*r/lam*NA).^6) - w;

plot(r,y)
