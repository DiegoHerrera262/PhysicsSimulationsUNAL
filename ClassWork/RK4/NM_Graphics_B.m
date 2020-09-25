%% Program for plotting Bessel function of task 2B
% Author: Diego Herrera
% Date: 14 - 09 - 20

%% Clear workspace
clear all; close all; clc;

%% Datafile name
InputDataFile = 'DataFiles/Task2B.txt';
OutputDataFile = 'Figures/Task2B.pdf';

%% Read data from file
SIR_Data = readtable(InputDataFile);
data = table2array(SIR_Data(:,1:2));

%% Start Plotting
myplot = figure(1);
% Plot susceptibles
plot(data(:,1),data(:,2),'LineWidth',2.0,'Color','b'); hold on;
%% Control aspect of graphics
xlabel('$\lambda$','FontSize',19,'interpreter','latex');
ylabel('$f(\lambda)$','FontSize',19,'interpreter','latex');
title('RK4 Integ. of Cost Function', '$f(\lambda) = R(r=1;\lambda)$',...
      'FontSize',19,'interpreter','latex');
%% Save figure
tag = input('Is the plot Proper to be saved? (1/0): ');
if tag == 1
  exportgraphics(myplot,OutputDataFile,'ContentType','vector');
end
