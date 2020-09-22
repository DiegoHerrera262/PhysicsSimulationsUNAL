%% Program that generates animation from dataset produced using C++ routines
% Author: Diego Alejandro Herrera Rojas
% Date: 20 - 09 - 20
% Description: This program takes as input some data file generated using
%              DiscreteElements tools programmed in C++. Produces a mpeg4
%              movie with simulated data. It works for 2D animations

function GenerateAnimation(AxisLayout,InputDataFile,OutputName)
  DOFs = 2;
  %% Read data from .txt file
  data = table2array(readtable(InputDataFile));
  dataSize = size(data);
  %% index for x coordinates
  idx = 2:DOFs:dataSize(2);
  idy = 3:DOFs:dataSize(2);
  %% Read time Values
  t = data(:,1);
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
  axis equal;
  axis(AxisLayout);
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
      movieFrames(int32(k/250)) = getframe(myframe);
    end
  end
  tag = input('Is the Animation proper to be saved? (1/0) : ');
  if tag == 1
    %% Save animations
    AnimationGenerator =  VideoWriter(OutputName,'MPEG-4');
    % Open writer, store video and close file
    open(AnimationGenerator);
    writeVideo(AnimationGenerator,movieFrames);
    close(AnimationGenerator);
  end
  clear mov
end
