function createGainGraph(X1, YMatrix1, Y1)
%CREATEFIGURE(X1, YMATRIX1, Y1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data
%  Y1:  vector of y data

%  Auto-generated by MATLAB on 04-Aug-2017 17:29:44

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'Parent',axes1);
set(plot1(1),'DisplayName','1010 nm');
set(plot1(2),'DisplayName','1020 nm');
set(plot1(3),'DisplayName','1030 nm');
set(plot1(4),'DisplayName','1050 nm');

% Create xlabel
xlabel('Pump Wavelength (nm)');

% Create ylabel
ylabel('Gain (m^{-1})');

% Create title
title('Gain and dQ/dt vs. Pump Wavelength');

%% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[980 1060]);
%% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[-50 50]);
grid(axes1,'on');
% Set the remaining axes properties
set(axes1,'YColor',[0 0 0],'YTick',[-33.33 -16.67 0 16.67 33.33 50]);
% Create axes
axes2 = axes('Parent',figure1,...
    'ColorOrder',[0.466 0.674 0.188;0.301 0.745 0.933;0.635 0.078 0.184;0 0.447 0.741;0.85 0.325 0.098;0.929 0.694 0.125;0.494 0.184 0.556]);
hold(axes2,'on');

% Create plot
plot(X1,Y1,'Parent',axes2,'DisplayName','dQ/dt');

% Create ylabel
ylabel('dQ/dt (mW/m)');

%% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes2,[980 1060]);
%% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes2,[-20 20]);
% Set the remaining axes properties
set(axes2,'Color','none','HitTest','off','YAxisLocation','right','YColor',...
    [0 0 0],'YTick',[-13.33 -6.67 0 6.67 13.33 20]);
% Create legend
legend(axes1,'show');
