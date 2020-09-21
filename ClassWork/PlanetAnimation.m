%% Program for making animations from file
% Author: Diego Alejandro Herrera
% Date: 20 - 09 - 20

%% Clear workspace
clear all; close all; clc;

%% Parameters for future function
DOFs = 2;
OutputName = 'BinarySystem';
InputDataFile = 'PrimeraSimulacion.txt';

%% Read data from .txt file
data = table2array(readtable(InputDataFile));
dataSize = size(data);
%% index for x coordinates
idx = 2:DOFs:dataSize(2);
idy = 3:DOFs:dataSize(2);
%% Read time Values
t = data(:,1);

%% Read Coordinates for Planets
% planet1_x = data(:,2); planet1_y = data(:,3);
% planet2_x = data(:,4); planet2_y = data(:,5);

%% Plot trajectories
myframe = figure(1);
set(gcf,'Position',[100,100,500,500]);
% Coordinates of animation elements
planets_x = data(1,idx);
planets_y = data(1,idy);
% Object that controls plot of data sets
dataPlots = scatter(planets_x,planets_y,500,'filled');
% Object that controls title of plot
titleString = ['Trajectories of Simulated Planets at t = ' num2str(t(1))];
animationTitle = title(titleString,'FontSize',19);
% Set axis properties
axis square;
axis([-2.5 2.5 -2.5 2.5]);
grid on;
% Axis labels
xlabel('X Coordinate','FontSize',19);
ylabel('Y Coordinate','Fontsize',19);

% Refresh plot
for k = 1:length(t)
  if mod(k,250) == 0
    planets_x = data(k,idx);
    planets_y = data(k,idy);
    set(dataPlots, 'XData', planets_x, 'YData', planets_y);
    titleString = ['Trajectories of Simulated Elements at t = ' num2str(t(k))];
    set(animationTitle,'String',titleString);
    % drawnow;
    movieFrames(int32(k/250)) = getframe(myframe);
    % Store frames for mp4 cration
    % movieFrames(k) = getframe;
  end
end

%% Save animations
AnimationGenerator =  VideoWriter(OutputName,'MPEG-4');
% Open writer, store video and close file
open(AnimationGenerator);
writeVideo(AnimationGenerator,movieFrames);
close(AnimationGenerator);

clear mov
