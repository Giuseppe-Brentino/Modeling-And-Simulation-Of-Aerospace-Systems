% Modeling and Simulation of Aerospace Systems
% Assignment #1
% Author: Giuseppe Brentino

%% Exercise 1

clearvars; close all; clc;
% Set up plot style 
plotStyle;

% Define symbolic variables and system of equations
x = sym('x', [1 2], 'real')';
f(x) = [x(2, :).^2 - x(1, :) - 2; -x(1, :).^2 + x(2, :) + 10];

% Find real solutions symbolically and cast to double
real_sol = solve(f);
real_sol = double(struct2array(real_sol))'; % [x_sol_1 , x_sol_2]

% Symbolically compute the Jacobian matrix
df = matlabFunction(jacobian(f, x), 'Vars', {x});

% Convert system of equations and Jacobian to MATLAB functions
f = matlabFunction(f, 'Vars', {x});

% Set tolerance and maximum number of iterations for Newton's method
tol = 1e-6;
n_max = 100;

% Initialize structures to store results
analytic = struct('x', zeros(2), 'max_e', zeros(2, 1), 'it', zeros(2, 1));
approx = analytic;
% Initialize initial guess matrix
x0 = zeros(2, 2);

% Iterate over each real solution
for i = 1:size(real_sol, 2)

    % Set initial guess for each solution
    x0(:, i) = real_sol(:, i) + ones(2, 1);

    % Analytic solution
    [analytic.x(:, i), analytic.max_e(i), analytic.it(i)] = newtonMethod(f, x0(:, i), tol, n_max, df);

    % Approximated solution 
    [approx.x(:, i), approx.max_e(i), approx.it(i)] = newtonMethod(f, x0(:, i), tol, n_max);

end

% Generate values for plotting
x1 = linspace(-5, 5, 100);
x2 = linspace(-5, 5, 100);
f1 = x2.^2 - 2;
f2 = x1.^2 - 10;

%%% Plot
figure()
hold on
grid on
% Plot zero contours of individual functions
plot(f1, x2, 'DisplayName', 'Zeros of $\textbf{f}_1$');
plot(x1, f2, 'DisplayName', 'Zeros of $\textbf{f}_2$');
% Plot found solutions and initial guesses
plot(analytic.x(1, :), analytic.x(2, :), 'o', 'DisplayName', 'Zeros of \textbf{f}', ...
    'MarkerFaceColor', "#EDB120")
plot(x0(1, :), x0(2, :), 'o', 'DisplayName', 'Initial guess', ...
    'MarkerFaceColor', "#7E2F8E")
% Labeling and legend
xlabel('$x_1$')
ylabel('$x_2$')
axis equal
xlim([-5, 10])
ylim([-10, 10])
legend('Orientation', 'vertical', 'Location', 'best');

%% Exercise 2
clearvars; close all; clc;

% Define the ODE and analytical solution
f = @(x, t) x - 2*t.^2 + 2;
x_an = @(t) 2*t.^2 + 4*t - exp(t) + 2;

% Define initial state and integration time
x0 = 1;
t0 = 0;
tf = 2;

% Set step sizes
h_vect = [0.5 0.2 0.05 0.01];
h_l = length(h_vect);

% Initialize cells and arrays to store results
T = cell(h_l, 2);
X = cell(h_l, 3);
cpu_time = zeros(n_iter, h_l, 2);
rel_error = cell(h_l, 2);

% Define number of iteration to have an accurate estimate of the cpu time
n_iter = 10000;

% Compute numerical solution and measure CPU time
for k = 1:n_iter
    for i = 1:h_l
        % RK2
        tic
        [T{i,1}, X{i,1}] = RK([t0, tf], h_vect(i), 2, f, x0);
        cpu_time(k, i, 1) = toc;
        % RK4
        tic
        [T{i,2}, X{i,2}] = RK([t0, tf], h_vect(i), 4, f, x0);
        cpu_time(k, i, 2) = toc;

        % Compare numerical and analytical solution (only on the first iteration)
        if k == 1
            X{i,3} = x_an(T{i,1})';
            rel_error{i,1} = abs(X{i,1} - X{i,3}) ./ X{i,3};
            rel_error{i,2} = abs(X{i,2} - X{i,3}) ./ X{i,3};
        end
    end
end

% Estimated CPU time
mean_cpu_time = zeros(h_l, 2);
for i = 1:h_l
    mean_cpu_time(i,:) = [mean(cpu_time(:, i, 1)); mean(cpu_time(:, i, 2))];
end

% Plots

% Relative error over time
figure()
grid on
for i = 1:h_l
    subplot(2, 2, i)
    semilogy(T{i,1}, rel_error{i,1}, 'DisplayName', 'RK2')
    hold on
    semilogy(T{i,1}, rel_error{i,2}, 'DisplayName', 'RK4')
    xlabel('Time [s]')
    ylabel('Relative error [-]')
    legend
    title("h = " + num2str(h_vect(i)));
end

% Time vs. maximum relative error
max_rel_err = zeros(4, 2);
for i = 1:h_l
    max_rel_err(i,1) = max(rel_error{i,1});
    max_rel_err(i,2) = max(rel_error{i,2}); 
end
figure()
grid on
semilogy(mean_cpu_time(:,1), max_rel_err(:,1), 'o', 'DisplayName', 'RK2')
hold on
semilogy(mean_cpu_time(:,2), max_rel_err(:,2), '*', 'DisplayName', 'RK4')
xlabel('CPU time [s]')
ylabel('Relative error [-]')
legend


%% functions
%%% general functions

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

%%% functions Ex. 1

function [x, err, i] = newtonMethod(f, x0, tol, n_max, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% newtonMethod solves a system of nonlinear equations using the Newton-Raphson
% method. It returns the solution x, the maximum error err, and the number 
% of iterations i.
%
%   Inputs:
%       - f: Function handle representing the system of equations.
%            Function handle, f(x), where x is a vector of real numbers.
%       - x0: Initial guess for the solution.
%             Column vector, size [n, 1], where n is the number of unknowns.
%       - tol: Tolerance for convergence.
%              Positive scalar.
%       - n_max: Maximum number of iterations.
%                Positive integer.
%       - varargin: Optional parameter for supplying the Jacobian matrix.
%                   Function handle, df(x), representing the Jacobian matrix.
%                   Size of df(x): [n, n], where n is the number of unknowns.
%
%   Outputs:
%       - x: Solution vector.
%            Column vector, size [n, 1], where n is the number of unknowns.
%       - err: Maximum error in the solution.
%              Positive scalar.
%       - i: Number of iterations performed.
%            Positive integer.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 0;
err = max(abs(f(x0)));

while err >= tol && i < n_max

    if nargin == 5
        % Use provided Jacobian matrix if available
        df = varargin{1}(x0);
    else
        % Numerical approximation of the Jacobian matrix using central
        % differences
        epsilon = max(sqrt(eps), sqrt(eps) * abs(x0));
        H = diag(epsilon);
        df = (2 * H) \ (f(x0 + H) - f(x0 - H));
    end

    % Update solution using Newton-Raphson method
    x = x0 - df \ f(x0);
    
    % Update initial guess, error and iteration counter
    x0 = x;
    err = max(abs(f(x)));
    i = i + 1;

end

% Check for convergence or reaching maximum iterations
if err >= tol || i == n_max
    error('The Method did not converge');
end

end

%%% functions Ex. 2

function [t_vect, x] = RK(T, h, method, f, x0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% RK solves an ordinary differential equation (ODE) using a
% general explicit Runge-Kutta method.
%
%   Inputs:
%       - T: Time span of the simulation.
%            Vector, [t_initial, t_final].
%       - h: Step size.
%            Positive scalar.
%       - method: Order of the Runge-Kutta method (2 or 4).
%            Scalar.
%       - f: Function handle representing the ODE.
%            Function handle, f(x, t), where x is the state vector and t is
%            the current time.
%       - x0: Initial state vector.
%            Column vector.
%
%   Outputs:
%       - t_vect: Time vector.
%                 Column vector.
%       - x: Solution matrix.
%            Each row represents the state vector at a specific time.
%

% define coefficients for RK2 (Heun) and RK4
    switch method
        case 2 
            alpha = [0; 1; 0];
            beta = [0, 0; 1, 0; 0.5, 0.5];
        case 4
            alpha = [0; 0.5; 0.5; 1; 1];
            beta = [0, 0, 0, 0;
                    0.5, 0, 0, 0; 
                    0, 0.5, 0, 0;
                    0, 0, 1, 0;
                    1/6, 1/3, 1/3, 1/6];
        otherwise 
            error('Choose a RK method between 2 and 4');
    end
    
    % initialize variables
    K = zeros(1, method);
    t_vect = T(1):h:T(end);
    steps = length(t_vect);
    x = [x0; zeros(steps-1, length(x0))];

    % Integrate
    for k = 1:steps-1

        % Compute each Runge-Kutta stage
        for i = 1:length(K)
            y = (x(k, :) + h * sum(beta(i, :) * K(:)));
            t = (t_vect(k) + alpha(i) * h);
            K(i) = f(y, t);
        end

        % Update the solution using the Runge-Kutta stages and coefficients
        s = 0;
        for i = 1:method
            s = s + beta(end, i) * K(i);
        end
        x(k+1, :) = x(k, :) + h * s;
        
    end
end








