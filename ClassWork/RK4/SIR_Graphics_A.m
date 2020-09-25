%% Programa para graficar resultados modelo SIR - RK4
% Author: Diego Herrera
% Date: 14 - 09 - 20

%% Clear workspace
clear all; close all; clc;

%% Datafile name
InputDataFile = 'DataFiles/Task1A.txt';
OutputDataFile = 'Figures/Task1A.pdf';

%% Read data from file
SIR_Data = readtable(InputDataFile);
data = table2array(SIR_Data(:,2:5));

%% Start Plotting
myplot = figure(1);
% Plot susceptibles
plot(data(:,1),data(:,2),'LineWidth',2.0,'Color','b'); hold on;
% Plot Infected
plot(data(:,1),data(:,3),'LineWidth',2.0,'Color','r'); hold on;
% Plot recovered
plot(data(:,1),data(:,4),'LineWidth',2.0,'Color','k');
%% Control aspect of graphics
xlabel('Time','FontSize',19,'interpreter','latex');
ylabel('Proportion of SIR','FontSize',19,'interpreter','latex');
title('Integration of SIR Model using RK4', '$\beta=0.35$ and $\gamma=0.08$',...
      'FontSize',19,'interpreter','latex');
%% Inlude text labels
text(55,0.2,'Inf.','FontSize',15,'Color','r','interpreter','latex');
text(21,0.9,'Susc.','FontSize',15,'Color','b','interpreter','latex');
text(40,0.85,'Recov.','FontSize',15,'Color','k','interpreter','latex');
%% Save figure
tag = input('Is the plot Proper to be saved? (1/0): ');
if tag == 1
  exportgraphics(myplot,OutputDataFile,'ContentType','vector');
end
