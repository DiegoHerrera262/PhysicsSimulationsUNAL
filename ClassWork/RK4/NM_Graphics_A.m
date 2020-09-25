%% Program for plotting Bessel function of task 2A
% Author: Diego Herrera
% Date: 14 - 09 - 20

%% Clear workspace
clear all; close all; clc;

%% Datafile name
InputDataFile = 'DataFiles/Task2A.txt';
OutputDataFile = 'Figures/Task2A.pdf';

%% Read data from file
SIR_Data = readtable(InputDataFile);
data = table2array(SIR_Data(:,2:4));

%% Start Plotting
myplot = figure(1);
% Plot susceptibles
plot(data(:,1),data(:,2),'LineWidth',2.0,'Color','b'); hold on;
%% Control aspect of graphics
xlabel('$r$','FontSize',19,'interpreter','latex');
ylabel('$R(r)$','FontSize',19,'interpreter','latex');
title('RK4 Solution of Drum Bessel', '$\lambda = 1$, $R(0) = 1$, $\frac{dR}{dr}(0) = 0$',...
      'FontSize',19,'interpreter','latex');
%% Save figure
tag = input('Is the plot Proper to be saved? (1/0): ');
if tag == 1
  exportgraphics(myplot,OutputDataFile,'ContentType','vector');
end
