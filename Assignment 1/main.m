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

% Define number of iteration to have an accurate estimate of the cpu time
n_iter = 10000;

% Initialize cells and arrays to store results
T = cell(h_l, 2);
X = cell(h_l, 3);
cpu_time = zeros(n_iter, h_l, 2);
rel_error = cell(h_l, 2);

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

%% Exercise 3

clearvars; close all; clc;
syms h
assume(h>=0);

% system RHS
A = @(alpha) [0 1; -1 2*cos(alpha)];

% operators F(h,alpha)
F_RK2 = @(h,alpha) eye( size(A(alpha),1) ) + h*A(alpha) + h^2/2*A(alpha)^2;
F_RK4 = @(h,alpha) h*A(alpha) + h^2/2*A(alpha)^2 + h^3/6*A(alpha)^3 + ...
    h^4/24*A(alpha)^4  + eye( size(A(alpha),1 ) );

% find largest timestep that guarantees stability
alpha = 0:0.01:pi;
h_RK2 = zeros(length(alpha),1);
h_RK4 = h_RK2;
lambda = h_RK2;
% alpha = pi
h_rk2_pi = double( solve( max( abs( eig( F_RK2(h,pi) ) ) )==1) );
h_RK2(end) = h_rk2_pi(h_rk2_pi~=0);
h_rk4_pi = double( solve( max( abs( eig( F_RK4(h,pi) ) ) )==1) );
h_RK4(end) = h_rk4_pi(h_rk4_pi~=0);
lambda(end) = max( eig( A(pi) ) );
% alpha  in [0, pi)
for i = length(alpha)-1:-1:1
    h_RK2(i) =  fzero(@(h)max( abs( eig( F_RK2(h,alpha(i)) ) ) )-1,h_RK2(i+1));
    h_RK4(i) =  fzero(@(h)max( abs( eig( F_RK4(h,alpha(i)) ) ) )-1,h_RK4(i+1));
    lambda(i) = max( eig( A(alpha(i)) ) );
end

figure
hold on
grid on
c=colororder('gem');
plot([real(h_RK2.*lambda) real(h_RK2.*lambda)], [imag(h_RK2.*lambda) -imag(h_RK2.*lambda)],...
    'Color',c(1,:),'DisplayName','RK2')
plot([real(h_RK4.*lambda) real(h_RK4.*lambda)], [imag(h_RK4.*lambda) -imag(h_RK4.*lambda)],...
    'Color',c(2,:),'DisplayName','RK4')
axis equal
xlabel('Re(h$\lambda$)')
ylabel('Im(h$\lambda$)')
legend('RK2','','RK4')

%% Ex 4
clearvars; close all; clc;

% problem data
x0 = [1;1];
t0 = 0;
tf = 1;
A = @(alpha) [0 1; -1 2*cos(alpha)];
tol_vect = [1e-3; 1e-4; 1e-5; 1e-6];
alpha_vect = linspace(0,pi,300);
t = [t0 tf];

% operators F(h,alpha)
F_RK1 = @(h,alpha) eye(size(A(alpha),1)) + h*(A(alpha));
F_RK2 = @(h,alpha) eye( size(A(alpha),1) ) + h*A(alpha) + h^2/2*A(alpha)^2;
F_RK4 = @(h,alpha) h*A(alpha) + h^2/2*A(alpha)^2 + h^3/6*A(alpha)^3 + ...
    h^4/24*A(alpha)^4  + eye( size(A(alpha),1 ) );

%init matrices
h_RK1 = zeros(length(alpha_vect),length(tol_vect));
h_RK2 = h_RK1;
h_RK4 = h_RK1;
fe_RK1 = zeros(length(tol_vect),1);
fe_RK2 = fe_RK1;
fe_RK4 = fe_RK1;
lambda = zeros(length(alpha_vect),1);

for i = 1:length(alpha_vect)
    %analitical solution
    x_an = expmv( A(alpha_vect(i)),x0,tf);
    % max eigenvalue
    lambda(i) = max(eig(A(alpha_vect(i))));

    for j = 1:length(tol_vect)

        opt = optimset('TolFun',1e-12);
        h_RK1(i,j)= fzero(@(h) objFcn(x_an, t, h, F_RK1, x0, tol_vect(j), alpha_vect(i)) ,[tol_vect(j)/2 tf], opt);
        h_RK2(i,j) = fzero(@(h) objFcn(x_an, t, h, F_RK2, x0, tol_vect(j), alpha_vect(i)) ,[tol_vect(j)/2 tf], opt);
        h_RK4(i,j) = fzero(@(h) objFcn(x_an, t, h, F_RK4, x0, tol_vect(j), alpha_vect(i)) ,[tol_vect(j)/2 tf], opt);

        % function evaluations vs alpha
        if i == length(alpha_vect)
            fun = @(x,t) A(alpha_vect(i))*x;
            [~, ~, fe_RK1(j)] = RK([t0, tf], h_RK1(i,j), 1, fun, x0);
            [~, ~, fe_RK2(j)] = RK([t0, tf], h_RK2(i,j), 2, fun, x0);
            [~, ~, fe_RK4(j)] = RK([t0, tf], h_RK4(i,j), 4, fun, x0);
        end

    end
end

%plot fe vs tol
figure
loglog(tol_vect, fe_RK1,'-o','DisplayName','RK1')
hold on
grid on
loglog(tol_vect, fe_RK2,'-o','DisplayName','RK2')
loglog(tol_vect, fe_RK4,'-o','DisplayName','RK4')
legend
xlabel('Tollerance')
ylabel('Function evaluations')

% plots h-lambda
figure %RK1
subplot(1,2,1)
hold on
grid on
plot([real(h_RK1(:,1).*lambda); nan; real(h_RK1(:,1).*lambda)], [imag(h_RK1(:,1).*lambda); nan; -imag(h_RK1(:,1).*lambda)] )
plot([real(h_RK1(:,2).*lambda); nan; real(h_RK1(:,2).*lambda)], [imag(h_RK1(:,2).*lambda); nan; -imag(h_RK1(:,2).*lambda)] )
xlabel('Re(h$\lambda$)')
ylabel('Im(h$\lambda$)')
axis equal
legend("tol = " + num2str(tol_vect(1:2)))
subplot(1,2,2)
hold on
grid on
plot([real(h_RK1(:,3).*lambda); nan; real(h_RK1(:,3).*lambda)], [imag(h_RK1(:,3).*lambda); nan; -imag(h_RK1(:,3).*lambda)] )
plot([real(h_RK1(:,4).*lambda); nan; real(h_RK1(:,4).*lambda)], [imag(h_RK1(:,4).*lambda); nan; -imag(h_RK1(:,4).*lambda)] )
xlabel('Re(h$\lambda$)')
ylabel('Im(h$\lambda$)')
axis equal
legend("tol = " + num2str(tol_vect(3:4)))

figure %RK2
hold on
grid on
plot([real(h_RK2(:,1).*lambda); nan; real(h_RK2(:,1).*lambda)], [imag(h_RK2(:,1).*lambda); nan; -imag(h_RK2(:,1).*lambda)] )
plot([real(h_RK2(:,2).*lambda); nan; real(h_RK2(:,2).*lambda)], [imag(h_RK2(:,2).*lambda); nan; -imag(h_RK2(:,2).*lambda)] )
plot([real(h_RK2(:,3).*lambda); nan; real(h_RK2(:,3).*lambda)], [imag(h_RK2(:,3).*lambda); nan; -imag(h_RK2(:,3).*lambda)] )
plot([real(h_RK2(:,4).*lambda); nan; real(h_RK2(:,4).*lambda)], [imag(h_RK2(:,4).*lambda); nan; -imag(h_RK2(:,4).*lambda)] )
xlabel('Re(h$\lambda$)')
ylabel('Im(h$\lambda$)')
axis equal
legend("tol = " + num2str(tol_vect))

figure %RK4
hold on
grid on
plot([real(h_RK4(:,1).*lambda); nan; real(h_RK4(:,1).*lambda)], [imag(h_RK4(:,1).*lambda); nan; -imag(h_RK4(:,1).*lambda)] )
plot([real(h_RK4(:,2).*lambda); nan; real(h_RK4(:,2).*lambda)], [imag(h_RK4(:,2).*lambda); nan; -imag(h_RK4(:,2).*lambda)] )
plot([real(h_RK4(:,3).*lambda); nan; real(h_RK4(:,3).*lambda)], [imag(h_RK4(:,3).*lambda); nan; -imag(h_RK4(:,3).*lambda)] )
plot([real(h_RK4(:,4).*lambda); nan; real(h_RK4(:,4).*lambda)], [imag(h_RK4(:,4).*lambda); nan; -imag(h_RK4(:,4).*lambda)] )
xlabel('Re(h$\lambda$)')
ylabel('Im(h$\lambda$)')
axis equal
legend("tol = " + num2str(tol_vect))

%% Ex 5
clearvars; close all; clc;

% Model matrix
A = @(alpha) [0 1; -1 2*cos(alpha)];

% BI2_theta operator
B_BI2 = @(alpha, theta, h) ( eye(size(A(alpha),1)) - (1-theta)*h*A(alpha) + 0.5*((1-theta)*h)^2 * A(alpha)^2) \...
    ( eye(size(A(alpha),1)) + theta*h*A(alpha) + 0.5*(theta*h)^2*A(alpha)^2 );

% define variables
syms h_sim
assume(h_sim>0.1);
theta = [0.1 0.3 0.4 0.7 0.9];
alpha = linspace(0,pi,300);
h = zeros(length(alpha), length(theta));
lambda = zeros(length(alpha),1);

% compute eigenvalues of A
for i = 1:length(alpha)
    lambda(i) = max( eig( A(alpha(i)) ) );
end

%plot stability regions
figure
c=colororder('gem');
axis equal
hold on
grid on

% compute the stability region for each theta
for j = 1:length(theta)

    %compute analitically the first timestep to have a good initial guess
    %for the next steps
    try
        alpha = linspace(0,pi,300);
        h(1,j) = double( solve( max( abs( eig( B_BI2(alpha(1),theta(j),h_sim) ) ) )==1) );
        flip_flag = false;
    catch
        alpha = linspace(pi,0,300);
        h(1,j) = double( solve( max( abs( eig( B_BI2(alpha(1),theta(j),h_sim) ) ) )==1) );
        flip_flag = true;
    end

    % compute h and lambda for each alpha between 0 and pi
    for i = 1:length(alpha)
        if j==1
            lambda(i) = max( eig( A(alpha(i)) ) );
        end
        if i ~= 1
            h0 = h(i-1,j);
            h(i,j) = fzero(@(h) max( abs( eig(B_BI2(alpha(i),theta(j),h)) ) )-1, h0 );
        end
    end
    % if necessary, flip the vector of h so to match the computed
    % eigenvalues lambda
    if flip_flag
        h(:,j) = flipud(h(:,j));
    end

    real_part = [real(h(:,j).*lambda);NaN; real(h(:,j).*lambda)];
    imag_part = [imag(h(:,j).*lambda);NaN; -imag(h(:,j).*lambda)];
    plot( real_part, imag_part,'DisplayName',"$BI2_{"+ num2str(theta(j)) + "}$",'Color',c(j,:));
end
legend
xlabel('Re(h$\lambda$)')
ylabel('Im(h$\lambda$)')

%%NOTA: PER THETA > 0.5 LA STABILITY REGION è DENTRO I CERCHI, PER THETA<0.5 E FUORI

%% Ex 6
clearvars; close all; clc;

%Problem data
B = [-180.5, 219.5; 179.5, -220.5];
x0 = [1; 1];
h = 0.1;
t0 = 0;
tf = 5;

% analytical solution
t = t0:h:tf;
x_an = expmv( B,x0,t);

% RK4
[~, x_rk4] = RK([t0, tf], h, 4, @(x,t)B*x, x0);

% plot AN-RK4
figure
subplot(1,2,1)
hold on
grid on
plot(t,x_an(1,:))
plot(t,x_rk4(:,1),'--')
xlabel('Time [s]')
ylabel('x1')
legend('Analytical','Numerical')
subplot(1,2,2)
hold on
grid on
plot(t,x_an(2,:))
plot(t,x_rk4(:,2),'--')
xlabel('Time [s]')
ylabel('x2')
legend('Analytical','Numerical')

%IEX4
[t,x_iex4] = IEX4 (@(x,t) B*x,h,[t0 tf],x0);

% plot AN-IEX4
figure
subplot(1,2,1)
hold on
grid on
plot(t,x_an(1,:))
plot(t,x_iex4(1,:),'--')
xlabel('Time [s]')
ylabel('x1')
legend('Analytical','Numerical')
subplot(1,2,2)
hold on
grid on
plot(t,x_an(2,:))
plot(t,x_iex4(2,:),'--')
xlabel('Time [s]')
ylabel('x2')
legend('Analytical','Numerical')

%iex4 error
err = abs( ( x_an - x_iex4 ) ./ x_an );
max(abs(X_AN-X_IEX))./max(abs(X_AN))
figure
semilogy(t,err(1,:))
hold on
grid on
xlabel('Time [s]')
ylabel('Normalized error')
%stability regions
A = @(alpha) [0 1; -1 2*cos(alpha)];
alpha = linspace(0,pi,300);
lambda = zeros(length(alpha),1);
h_rk4 = lambda;
h_iex4 = lambda;

% RK4 operator
F_RK4 = @(h,alpha) h*A(alpha) + h^2/2*A(alpha)^2 + h^3/6*A(alpha)^3 + ...
    h^4/24*A(alpha)^4  + eye( size(A(alpha),1 ) );
%IEX4 operator
F_IEX4=@(h,alpha) -(1/6)*(eye(size(A(alpha),1)) -(h)*A(alpha))^(-1)...
    +4*(eye(size(A(alpha),1)) -(h/2)*A(alpha))^(-2) ...
    -(27/2)*(eye(size(A(alpha),1)) -(h/3)*A(alpha))^(-3)...
    +(32/3)*(eye(size(A(alpha),1)) -(h/4)*A(alpha))^(-4);

syms h_sim real
assume(h_sim>0)
h_rk4(end) =  double( solve( max( abs( eig( F_RK4(h_sim,pi) ) ) )==1));
h_iex4(1) = max(double( solve( max( abs( eig( F_IEX4(h_sim,0) ) ) )==1)));

for i = 1:length(alpha)
    lambda(i) = max( eig( A(alpha(i)) ) );
end

for i = 2:length(alpha)
    h_iex4(i) = fzero(@(h) max(abs(eig(F_IEX4(h,alpha(i)))))-1,h_iex4(i-1));
end
for i = length(alpha)-1:-1:1
    h_rk4(i) = fzero(@(h) max(abs(eig(F_RK4(h,alpha(i)))))-1,h_rk4(i+1));
end

figure
c=colororder('gem');
hold on
%IEX4
iex4_stability_r = real(lambda.*h_iex4);
iex4_stability_i = imag(lambda.*h_iex4);
patch([-42 -42 15 15],[-10 10 10 -10],c(1,:),'FaceAlpha',0.2,'EdgeAlpha',0, ...
    'DisplayName','IEX4 stability');
patch(iex4_stability_r, iex4_stability_i, [1 1 1],'EdgeColor',c(1,:),...
    'EdgeAlpha',0,'FaceAlpha',1,'HandleVisibility','off')
patch(iex4_stability_r, -iex4_stability_i,[1 1 1],'EdgeColor',c(1,:),...
    'EdgeAlpha',0,'FaceAlpha',1,'HandleVisibility','off')
plot([real(lambda.*h_iex4); nan; real(lambda.*h_iex4)],...
    [imag(lambda.*h_iex4); nan; -imag(lambda.*h_iex4)],'DisplayName','IEX4')

%RK4
rk4_stability_r = real(lambda.*h_rk4);
rk4_stability_i = imag(lambda.*h_rk4);
patch(rk4_stability_r, rk4_stability_i, c(2,:),'EdgeColor',c(2,:),...
    'EdgeAlpha',0,'FaceAlpha',0.2,'DisplayName','RK4 stability region')
patch(rk4_stability_r, -rk4_stability_i, c(2,:),'EdgeColor',c(2,:),...
    'EdgeAlpha',0,'FaceAlpha',0.2,'HandleVisibility','off')
plot([rk4_stability_r; nan; rk4_stability_r], [rk4_stability_i; nan; -rk4_stability_i],'DisplayName','RK4')
%Problem eigs
lambda_prob = h*eig(B);
plot(real(lambda_prob),imag(lambda_prob),'.','MarkerSize',25,'DisplayName','h$\lambda$ of the problem')
axis equal
grid off
xlim([-40 14]);
ylim([-8 8])

legend

%% Ex 7
clearvars; close all; clc;

%data
h = 0.1;
t0 = 0;
tf = 3;
x0 = [1;1];


% compute eigenvalues of the problem
syms t real
assume(t>=0);
syms x(t)
% compute x1 analytically
x1_sym = dsolve(diff(x,1,t)==-(5/2)*(1+8*sin(t))*x(t),x(0)==1);

% Linear formulation of the problem for each t
F = [-(5/2)*(1+8*sin(t)), 0;
    1, (1-x1_sym)];

time = t0:h:tf;
lambda_fun = matlabFunction( eig(F) );
lambda = lambda_fun(time);

figure
hold on
grid on
plot(time,lambda(1,:), time,lambda(2,:))
xlabel('Time [s]')
ylabel('$\lambda$')
legend('$\lambda_1$', '$\lambda_2$')

figure
c=colororder('gem');
hold on
grid on
plot(time,h*lambda(1,:),'DisplayName','$h\lambda_1$')
plot(time,h*lambda(2,:),'DisplayName','$h\lambda_2$')
xlabel('Time [s]')
ylabel('$h\lambda$')
yline(-0.5495,'--','LineWidth',1.5,'Color',c(3,:),'DisplayName','AB3 min $h\lambda$')
yline(-6,'--','LineWidth',1.5,'Color',c(4,:),'DisplayName','AM3 min $h\lambda$')
yline(-1.7327,'--','LineWidth',1.5,'Color',c(5,:),'DisplayName','ABM3 min $h\lambda$')
yline(6.6434,'--','LineWidth',1.5,'Color',c(6,:),'DisplayName','BDF3 max $h\lambda$')
ylim([-7 7])
xlim([0 3])
legend


% AB3
[time,x_AB3] = AB3(@f_ex7,[t0,tf],x0,h);
figure
subplot(1,2,1)
plot(time,x_AB3(1,:))
subplot(1,2,2)
plot(time,x_AB3(2,:))

%AM3
[time,x_AM3] = AM3(@f_ex7,[t0,tf],x0,h);
figure
subplot(1,2,1)
plot(time,x_AM3(1,:))
subplot(1,2,2)
plot(time,x_AM3(2,:))


%ABM3
[time,x_ABM3] = ABM3(@f_ex7,[t0,tf],x0,h);
figure
subplot(1,2,1)
plot(time,x_ABM3(1,:))
subplot(1,2,2)
plot(time,x_ABM3(2,:))

%BDF3

[time,x_BDF3] = BDF3(@f_ex7,[t0,tf],x0,h);
figure
subplot(1,2,1)
plot(time,x_BDF3(1,:))
subplot(1,2,2)
plot(time,x_BDF3(2,:))

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

function [t_vect, x, fe] = RK(T, h, method, f, x0)
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
%       - fe: Number of function evaluations. Scalar
%
% define coefficients for RK1(Forward Euler), RK2 (Heun) and RK4

arguments
    T (2,1)
    h (1,1)
    method (1,1)
    f
    x0 (1,:)
end

switch method
    case 1
        alpha = [0; 1];
        beta = [0;1];
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
        error('Choose a RK method between 1, 2 and 4');
end

% initialize variables
K = zeros(length(x0), method);
t_vect = T(1):h:T(end);
steps = length(t_vect);
x = [x0; zeros(steps-1, length(x0))];
fe = 0;

% Integrate
for k = 1:steps-1

    % Compute each Runge-Kutta stage
    for i = 1:size(K,2)
        if i~=1
            y = x(k, :) + h * sum(beta(i, :) .* K(:,i-1),2)';
        else
            y = x(k, :);
        end
        t = (t_vect(k) + alpha(i) * h);
        K(:,i) = f(y', t);
        fe = fe+1;

    end


    % Update the solution using the Runge-Kutta stages and coefficients
    s = zeros(1,length(x0));
    for i = 1:method
        s = s + (beta(end, i) * K(:,i))';
    end
    x(k+1, :) = x(k, :) + h * s;

end
end

%%% function Ex. 1

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

%%% function Ex. 4
function J = objFcn(x_an, t, h, f, x0, tol, alpha)
f = f(h,alpha);
n_steps = (t(2) - t(1))/h ;
J =  max( abs(x_an-f^(n_steps)*x0) ) - tol;
end

%%% function Ex. 6

function [t,x] = IEX4 (f,h,t,x0)

alpha = [-1/6 4 -27/2 32/3];
opt = optimset('Display','off');

t = t(1):h:t(2);
x = zeros(length(x0),length(t));
x(:,1) = x0;

for i = 2:length(t)
    %1st predictor
    k1 = fsolve(@(k) h*f(k, t(i)+h) - k  + x(:,i-1), x(:,i-1), opt );

    %2nd predictor
    k2a = fsolve(@(k) h/2*f(k, t(i-1)+h/2) - k  + x(:,i-1), x(:,i-1), opt );
    k2 = fsolve(@(k) h/2*f(k, t(i-1)+h/2) - k  + k2a, x(:,i-1), opt );

    %3rd predictor
    k3a = fsolve(@(k) h/3*f(k, t(i-1)+h/3) - k  + x(:,i-1), x(:,i-1), opt );
    k3b = fsolve(@(k) h/3*f(k, t(i-1)+h/3) - k  + k3a, x(:,i-1), opt );
    k3 = fsolve(@(k) h/3*f(k, t(i-1)+h/3) - k  + k3b, x(:,i-1), opt );

    %4th predictor
    k4a = fsolve(@(k) h/4*f(k, t(i-1)+h/4) - k  + x(:,i-1), x(:,i-1), opt );
    k4b = fsolve(@(k) h/4*f(k, t(i-1)+h/4) - k  + k4a, x(:,i-1), opt );
    k4c = fsolve(@(k) h/4*f(k, t(i-1)+h/4) - k  + k4b, x(:,i-1), opt );
    k4 = fsolve(@(k)  h/4*f(k, t(i-1)+h/4) - k  + k4c, x(:,i-1), opt );

    % final sol
    x(:,i) = alpha(1)*k1 + alpha(2)*k2 + alpha(3)*k3 + alpha(4)*k4;
end

end

%%% functions Ex. 7

function dx = f_ex7 (x,t)

dx = [-5/2 * ( 1+8*sin(t) ) * x(1);
    ( 1-x(1) )*x(2) + x(1)];
end

function [t,x] = AB3(f,time,x0,h)
t0 = time(1);
tf = time(2);
t = t0:h:tf;
x = zeros(length(x0),length(time));

x(:,1) = x0;
for k = 1:length(t)-1
    %use RK4 for startup
    if k < 3
        [~, x_t, ~] = RK([t(k), t(k)+h], h, 4, f, x(:,k));
        x(:,k+1) = x_t(end,:)';
    else
        x(:,k+1) = x(:,k) + h/12 * ( 23*f(x(:,k),t(k)) -16*f(x(:,k-1),t(k-1))...
            +5*f(x(:,k-2),t(k-2)) );
    end
end

end

function [t,x] = AM3(f,time,x0,h)
t0 = time(1);
tf = time(2);
t = t0:h:tf;
x = zeros(length(x0),length(time));
opt = optimset('Display','off');
x(:,1) = x0;
for k = 1:length(t)-1
    %use RK4 for startup
    if k < 2
        [~, x_t, ~] = RK([t(k), t(k)+h], h, 4, f, x(:,k));
        x(:,k+1) = x_t(end,:)';
    else
        fun = @(y) x(:,k) + h/12 * ( 5*f(y,t(k+1)) +8*f(x(:,k),t(k)) -f(x(:,k-1),t(k-1)) ) -y;
        x(:,k+1) = fsolve(fun,x(:,k),opt);
    end
end

end

function [t,x] = ABM3(f,time,x0,h)
t0 = time(1);
tf = time(2);
t = t0:h:tf;
x = zeros(length(x0),length(time));
x(:,1) = x0;

for k = 1:length(t)-1
    % startup with rk4
    if k < 3
        [~, x_t, ~] = RK([t(k), t(k)+h], h, 4, f, x(:,k));
        x(:,k+1) = x_t(end,:)';
    else
        % Predictor (AB3)
        x(:,k+1) = x(:,k) + h/12 * ( 23*f(x(:,k),t(k)) -16*f(x(:,k-1),t(k-1))...
            +5*f(x(:,k-2),t(k-2)) );
        % Corrector (AM3)
        x(:,k+1) = x(:,k) + h/12*( 5*f(x(:,k+1),t(k+1)) +8*f(x(:,k),t(k)) - f(x(:,k-1),t(k-1)) );
    end

end

end


function [t,x] = BDF3(f,time,x0,h)
t0 = time(1);
tf = time(2);
t = t0:h:tf;
x = zeros(length(x0),length(time));
opt = optimset('Display','off');
x(:,1) = x0;
for k = 1:length(t)-1
    %use RK4 for startup
    if k < 3
        [~, x_t, ~] = RK([t(k), t(k)+h], h, 4, f, x(:,k));
        x(:,k+1) = x_t(end,:)';
    else
        fun = @(y) 18/11*x(:,k) -9/11*x(:,k-1) +2/11*x(:,k-2) + 6/11*h*f(y,t(k+1)) - y;
        x(:,k+1) = fsolve(fun,x(:,k),opt);
    end
end
end
