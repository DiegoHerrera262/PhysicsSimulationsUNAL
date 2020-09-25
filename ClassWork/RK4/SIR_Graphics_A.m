%% Programa para graficar resultados modelo SIR - RK4
% Author: Diego Herrera
% Date: 14 - 09 - 20

%% Clear workspace
clear all; close all; clc;

%% Datafile name
InputDataFile = 'DataFiles/Punto1A.txt';
OutputDataFile = 'Figures/Punto1A.pdf';

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
xlabel('Time','FontSize',19);
ylabel('Proportion of SIR','FontSize',19);
title('Integration of SIR Model using RK4','FontSize',19);
%% Save figure
exportgraphics(myplot,OutputDataFile,'ContentType','vector');
