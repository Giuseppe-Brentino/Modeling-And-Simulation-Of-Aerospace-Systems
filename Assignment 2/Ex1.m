% Modeling and Simulation of Aerospace Systems (2023/2024)
% Assignment #2: Exercise 1
% Author: Giuseppe Brentino

close all; clearvars; clc;
%% Part1 : Thermal properties

properties.A = 1; %m^2 unitary area
% Tantalum carbide
properties.k1 = 40;   % W/(m K) volendo si può interpolare curva
properties.l1 = 5e-4; %m
properties.R1 = properties.l1/properties.k1;

% Graphite
properties.c2 = 8110; %J/(kg*K) (NIST)
properties.k2 = 75;      % W/(m K) si può interpolare curva
properties.rho2 = 2100; %g/cm^3 wikipedia
properties.l2 = 0.01; 
properties.R2 = properties.l2/properties.k2;
properties.C2 = properties.rho2*properties.l2*properties.c2;

%Phenolic resin SPGB0
properties.rho4 = 1340; %kg/m^3 
properties.c4 = 1250; % J/(kg K);
properties.k4 = 0.26; % W/(m K)
properties.l4 = 0.002; 
properties.R4 = properties.l4/properties.k4;
properties.C4 = properties.rho4*properties.l4*properties.c4;

% interface
properties.R3 = (properties.R2 + properties.R4)/2; 

% Aluminum
properties.k5 = 237; %W/(m K); https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10144406/#:~:text=Aluminum%20has%20a%20thermal%20conductivity,1%2C2%2C3%5D.
properties.l5 = 0.001;
properties.R5 = properties.l5/properties.k5;

%interface
properties.Ti = @(t) 20 + (980/1.*t).*(t<=1) + 980.*(t>1) + 273.15;
properties.To = 293.15;

%% Part 1.2 causal modeling

x0 = 293.15*ones(5,1);
ode_tol = 1e-6;
opt = odeset('RelTol',ode_tol,'AbsTol',ode_tol);
[t,T] = ode15s(@thermalModel,[0 60],x0,opt,properties);

figure()
hold on
grid on
plot(t,T(:,1))
plot(t,T(:,2))
plot(t,T(:,3))
plot(t,T(:,4))
plot(t,T(:,5))
legend('T1','T2','T3','T4','T5')

%% Part 1.3: acasual modeling

% model = sim('modello.slx');


%% functions

function dx = thermalModel(t,x,properties)
To = properties.To;
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
