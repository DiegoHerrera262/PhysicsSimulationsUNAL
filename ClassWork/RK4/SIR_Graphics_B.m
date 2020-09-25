%% Programa para graficar Modelo SIR exacto
% Author: Diego Herrera
% Date: 14 - 09 - 20

%% Limpiar workspace
clear all; close all; clc;

OutputDataFile = 'Figures/Task1B.pdf';

%% Definicion de Funciones Importantes
% Parametros del modelo
R0 = 0.0; S0 = 0.999; I0 = 0.001;
gamma = 0.08; beta = 0.35; rho = gamma/beta;
% Variables de solucion analitica
u = logspace(log(0.155),log(1.0),1500);  % variable u de la solucion analitica
S = @(x) S0 * x;                        % Definicion de Susceptibles
R = @(x) R0 - rho * log(x);             % Definicion de Recuperados
I = @(x) 1 - S(x) - R(x);               % Definicion de integral Infectados
tinteg = @(x) 1.0 ./ beta .* 1.0 ./(x .* I(x));
t = @(x) integral(tinteg,x,1);  % Definicion de infectados

%% Calculo variables modelo
Sdata = S(u);
Rdata = R(u);
Idata = I(u);
tdata = zeros(1,length(u));
for idx = 1:length(u)
    tdata(idx) = t(u(idx));
end

%% Realiza graficas
myplot = figure(1);
inrange = (tdata >= 0) & (tdata <= 80);
plot(tdata(inrange),Sdata(inrange),'LineWidth',2.0,'Color','b'); hold on;
plot(tdata(inrange),Rdata(inrange),'LineWidth',2.0,'Color','k'); hold on;
plot(tdata(inrange),Idata(inrange),'LineWidth',2.0,'Color','r'); hold on;
%% Control aspect of graphics
xlabel('Time','FontSize',19,'interpreter','latex');
ylabel('Proportion of SIR','FontSize',19,'interpreter','latex');
title('Analytical Solution of SIR Model', '$\beta=0.35$ and $\gamma=0.08$',...
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
