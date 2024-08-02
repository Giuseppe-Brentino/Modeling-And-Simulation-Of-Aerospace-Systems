clearvars; close all; clc;
plotStyle;
%% Part 1
load EX_2_1_10s.mat

% angular speed
figure
hold on
grid on
plot(data(:,1),data(:,2))
xlabel('Time [s]')
ylabel('$\omega$ [rad/s]')

% Voltage
figure
hold on
grid on
plot(data(:,1),data(:,3))
xlabel('Time [s]')
ylabel('Voltage [V]')

load EX_2_1_120s.mat

% Temperature
figure
hold on
grid on
plot(data(:,1),data(:,2))
xlabel('Time [s]')
ylabel('Temperature [K]')

%heat flux
figure
hold on
grid on
plot(data(:,1),data(:,3))
xlabel('Time [s]')
ylabel('Heat flux [W]')

%% Part 2
load EX2_2.mat

% Gearbox Temperature
figure
hold on
grid on
plot(data(:,1),data(:,3),'DisplayName','Gearbox temperature')
xlabel('Time [s]')
ylabel('Temperature [°C]')
yline(60,'k--','LineWidth',1.5,'DisplayName','Max T')
yline(40,'r--','LineWidth',1.5,'DisplayName','Min T')
legend

%Tank level
figure
hold on
grid on
xlabel('Time [s]')
yyaxis left
plot(data(:,1),data(:,2),'DisplayName','Water level')
ylabel('Water level [m]')
yyaxis right
plot(data(:,1),data(:,6)*3600,'DisplayName','Water flow')
ylabel('Water flow [$m^3/h$]')
legend

% Sink tank temperature
figure
hold on
grid on
xlabel('Time [s]')
yyaxis left
plot(data(:,1),data(:,5)-273.15,'DisplayName','Water temperature')
ylabel('Temperature [°C]')
yline(10,'k--','LineWidth',1.5,'DisplayName','Max T')
yline(5,'r--','LineWidth',1.5,'DisplayName','Min T')
yyaxis right
plot(data(:,1),data(:,4),'DisplayName','Heat flux')
ylabel('Heat Flux [W]')
legend

%% function
function plotStyle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set figure properties for better looking plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interpreter:
set(0, 'defaultTextInterpreter', 'Latex')
set(0, 'defaultAxesTickLabelInterpreter', 'Latex')
set(0, 'defaultLegendInterpreter', 'Latex')
set(0,'defaultAxesXGrid','on')
set(0,'defaultAxesYGrid','on')
% lines:
set(0,'defaultLineLineWidth', 1.5);
set(0,'defaultLineMarkerSize',6) ;
% legend:
set(0, 'defaultLegendLocation','southoutside');
set(0, 'defaultLegendOrientation','horizontal');
set(0, 'defaultLegendFontSize',12);
% axes:
set(0,'defaultAxesFontSize',16);
end
