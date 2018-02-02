clear all

delT = 0:0.1:200;
hRAW = xlsread('h_values.xlsx');
[sorted_h, sortIndex_h] = unique(hRAW(:,1));
hRAW = [sorted_h hRAW(sortIndex_h,2)];
%plot(hRAW(:,1),hRAW(:,2))
hINTERP = [delT; interp1(hRAW(:,1),hRAW(:,2),delT)].';
plot(delT,hINTERP(:,2))

csvwrite('h_values.csv',hINTERP);