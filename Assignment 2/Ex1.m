% Modeling and Simulation of Aerospace Systems (2023/2024)
% Assignment #2: Exercise 1
% Author: Giuseppe Brentino

close all; clearvars; clc;
%% Part1 : Thermal properties

properties.A = 1; %m^2 unitary area
% Tantalum carbide
properties.k1 = 40;   % W/(m K)
properties.l1 = 5e-4; %m
properties.R1 = properties.l1/properties.k1;

% Graphite
properties.c2 = 8110; %J/(kg*K) (NIST)
properties.k2 = 75;      % W/(m K) 
properties.rho2 = 2100; % Kg/m^3
properties.l2 = 0.02; 
properties.R2 = properties.l2/properties.k2;
properties.C2 = properties.rho2*properties.l2*properties.c2;


%Phenolic resin SPGB0
properties.rho4 = 1340; %kg/m^3 
properties.c4 = 1250; % J/(kg K);
properties.k4 = 0.26; % W/(m K)
properties.l4 = 0.005; 
properties.R4 = properties.l4/properties.k4;
properties.C4 = properties.rho4*properties.l4*properties.c4;

% interface
properties.R3 = (properties.R2 + properties.R4)/2; 

% Aluminum
properties.k5 = 237; %W/(m K); 
properties.l5 = 0.001;
properties.R5 = properties.l5/properties.k5;

% system properties
properties.Tgas = 1273.15; %[K]
properties.T0 = 293.15;
properties.Ti = @(t) properties.T0 + ((properties.Tgas-properties.T0)*t).*(t<=1)...
    + (properties.Tgas-properties.T0).*(t>1);


%% Part 1.3: acasual modeling

x0 = 293.15*ones(5,1);
ode_tol = 1e-11;
tf = 60;

model = sim('model.slx');

% plot one-node-model temperature
figure
hold on 
grid on
plot(model.tout,model.Ti,'k--')
plot(model.tout,model.one_node.T_layers.T1.Data(:,1))
plot(model.tout,model.one_node.T_layers.T2.Data(:,1))
plot(model.tout,model.one_node.T_layers.T3.Data(:,1))
plot(model.tout,model.one_node.T_layers.T4.Data(:,1))
plot(model.tout,model.one_node.T_layers.T5.Data(:,1))
legend('$T_i$','$T_1$','$T_2$','$T_3$','$T_4$','$T_5$')
xlabel('Time [s]')
ylabel('Temperature [K]')

% plot two-node-model temperature
figure
hold on 
grid on
plot(model.tout,model.Ti,'k--')
plot(model.tout,model.two_nodes.T_layers.T_1.Data(:,1))
plot(model.tout,model.two_nodes.T_layers.T_21.Data(:,1))
plot(model.tout,model.two_nodes.T_layers.T_22.Data(:,1))
plot(model.tout,model.two_nodes.T_layers.T_3.Data(:,1))
plot(model.tout,model.two_nodes.T_layers.T_41.Data(:,1))
plot(model.tout,model.two_nodes.T_layers.T_42.Data(:,1))
plot(model.tout,model.two_nodes.T_layers.T_5.Data(:,1))
legend('$T_i$','$T_1$','$T_{21}$','$T_{22}$','$T_3$','$T_{41}$','$T_{42}$','$T_5$')
xlabel('Time [s]')
ylabel('Temperature [K]')

% Two-nodes vs 1-node models
tf = 600;
modelLong = sim('model.slx');
deltaT1 = modelLong.two_nodes.T_interfaces.T_int12.Data(:,1) - modelLong.one_node.T_interfaces.T_int12.Data(:,1);
deltaT2 = modelLong.two_nodes.T_interfaces.T_int23.Data(:,1) - modelLong.one_node.T_interfaces.T_int23.Data(:,1);
deltaT3 = modelLong.two_nodes.T_interfaces.T_int34.Data(:,1) - modelLong.one_node.T_interfaces.T_int34.Data(:,1);
deltaT4 = modelLong.two_nodes.T_interfaces.T_int45.Data(:,1) - modelLong.one_node.T_interfaces.T_int45.Data(:,1);

figure
hold on
grid on
plot(modelLong.tout,deltaT1)
plot(modelLong.tout,deltaT2)
plot(modelLong.tout,deltaT3)
plot(modelLong.tout,deltaT4)
legend('$T_{1-2}$','$T_{2-3}$','$T_{3-4}$','$T_{4-5}$')
xlabel('Time [s]')
ylabel('$T_{2nodes}-T_{1node}$ [K]')
%% Part 1.2 causal modeling

opt = odeset('RelTol',ode_tol,'AbsTol',ode_tol);

[t,T] = ode15s(@thermalModel,model.tout,x0,opt,properties);


% casual modeling plot
figure()
hold on
grid on
plot(t,properties.Ti(t),'k--');
plot(t,T(:,1))
plot(t,T(:,2))
plot(t,T(:,3))
plot(t,T(:,4))
plot(t,T(:,5))
legend('$T_i$','$T_1$','$T_2$','$T_3$','$T_4$','$T_5$') 
xlabel('Time [s]')
ylabel('Temperature [K]')

% plot acasual vs causal modeling
error.T1 = abs( model.one_node.T_layers.T1.Data(:,1) - T(:,1) ) ./ T(:,1);
error.T2 = abs( model.one_node.T_layers.T2.Data(:,1) - T(:,2) ) ./ T(:,2);
error.T3 = abs( model.one_node.T_layers.T3.Data(:,1) - T(:,3) ) ./ T(:,3);
error.T4 = abs( model.one_node.T_layers.T4.Data(:,1) - T(:,4) ) ./ T(:,4);
error.T5 = abs( model.one_node.T_layers.T5.Data(:,1) - T(:,5) ) ./ T(:,5);

figure
semilogy(t,error.T1);
hold on
grid on
semilogy(t,error.T2);
semilogy(t,error.T3);
semilogy(t,error.T4);
semilogy(t,error.T5);
xlabel('Time [s]')
ylabel('Normalized relative error')
legend('$T_1$','$T_2$','$T_3$','$T_4$','$T_5$')
%% functions

function dx = thermalModel(t,x,properties)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% thermalModel - Differential thermal model of the system.
%
% Inputs:
%   t          - Current time
%   x          - State vector [T1, T2, T3, T4, T5]
%   properties - Structure containing the thermal properties:
%                .T0 - Ambient temperature
%                .Ti - Function handle for the input temperature (Ti(t))
%                .R1, R2, R3, R4, R5 - Thermal resistances
%                .C2, C4 - Thermal capacitances
%
% Outputs:
%   dx         - Time derivative of the state vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

To = properties.T0;
Ti = properties.Ti(t);

R1 = properties.R1;
R2 = properties.R2;
R3 = properties.R3;
R4 = properties.R4;
R5 = properties.R5;

C2 = properties.C2;
C4 = properties.C4;

dx(1) = ( Ti-x(1) )/(0.5*R1) - ( x(1)-x(2) )/( 0.5*(R1+R2) );
dx(2) = 1/C2 * ( (x(1)-x(2))/(0.5*(R1+R2)) - (x(2)-x(3))/(0.5*(R2+R3)) );
dx(3) = (x(2)-x(3))/(0.5*(R2+R3)) - (x(3)-x(4))/(0.5*(R3+R4));
dx(4) = 1/C4 * ( (x(3)-x(4))/(0.5*(R3+R4)) - (x(4)-x(5))/(0.5*(R4+R5)) );
dx(5) = (x(4)-x(5))/(0.5*(R4+R5)) - (x(5)-To)/(0.5*R5);

dx = dx';
end

