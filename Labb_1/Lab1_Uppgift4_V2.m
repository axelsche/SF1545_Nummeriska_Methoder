%%% The code models the trend in USDSEK Forex chart over the past 1000+
%%% days. Ititialy with a linear fit to later use a model and finaly solve
%%% it a third way using gauss-Newtons method. From visual annalysis the
%%% period appers to  be 919 days, but use 920 in the guess, whih later
%%% is shown to be 1143 from G-N method


% Close all figures, clear all variables from workspace, and clear command window
close all;
clear variables;
clc;

%% Load Data
% Load the data from 'dollarkurs.mat', which should contain variables 'day' and 'USDSEK'
load('dollarkurs.mat');  % Ensure 'day' and 'USDSEK' variables are in the file

% Get the number of data points
N = length(day);  % Number of data points

%% Task a: Linear Model Fitting
% Task a: Fit a linear model to the data: USDSEK = c0 + c1 * day
% Construct the design matrix A for linear regression with intercept (c0) and slope (c1)
A = [ones(N, 1), day];

% Solve the normal equations to find the least squares solution for coefficients c0 and c1
coeffs = A \ USDSEK;
c0 = coeffs(1);  % Intercept
c1 = coeffs(2);  % Slope

% Compute the fitted values using the linear model
g_linear = A * coeffs;

% Compute the mean squared error (MSE) of the linear model
E_linear = mean((USDSEK - g_linear).^2);

% Display the results
fprintf('TASK A - Linear Model Coefficients:\n');
fprintf('c0 = %.6f, c1 = %.6f\n', c0, c1);
fprintf('Mean Squared Error: E = %.6f\n\n', E_linear);

% Plot the original data and the fitted linear model
figure;
plot(day, USDSEK, 'bo', 'DisplayName', 'Data');  % Original data points
hold on;
plot(day, g_linear, 'r-', 'LineWidth', 2, 'DisplayName', 'Linear Model');  % Linear fit
xlabel('Days');
ylabel('USD/SEK');
title('Linear Model Fitting to Dollar Exchange Rate');
legend('Location', 'best');
grid on;
hold off;

%% Task b: Linear + Periodic Model Fitting
% Task b: Fit a model that includes both linear and periodic components
% The model is: USDSEK = d0 + d1 * day + d2 * sin(2*pi*day/L) + d3 * cos(2*pi*day/L)
% where L is the period of the periodic component
% We will use a fixed period L for this task

% Compute the residuals from the linear model
residuals_linear = USDSEK - g_linear;

% Plot the residuals to observe any periodic behavior
figure;
plot(day, residuals_linear, 'b-');
xlabel('Days');
ylabel('Residuals (USD/SEK)');
title('Residuals Between Linear Model and Data');
grid on;

% Initial guess for the period L of the periodic component (e.g., from observing the residual plot)
L = 920;  % Initial guess for the period (in days)

% Construct the design matrix A_periodic including the periodic components
A_periodic = [ones(N, 1), day, sin(2 * pi * day / L), cos(2 * pi * day / L)];

% Solve the normal equations to find the least squares solution for coefficients d0, d1, d2, d3
coeffs_periodic = A_periodic \ USDSEK;
d0 = coeffs_periodic(1);  % Intercept
d1 = coeffs_periodic(2);  % Linear term coefficient
d2 = coeffs_periodic(3);  % Coefficient for sine term
d3 = coeffs_periodic(4);  % Coefficient for cosine term

% Compute the fitted values using the linear + periodic model
g_periodic = A_periodic * coeffs_periodic;

% Compute the mean squared error (MSE) of the linear + periodic model
E_periodic = mean((USDSEK - g_periodic).^2);

% Display the results
fprintf('TAKS B - Linear + Periodic Model Coefficients:\n');
fprintf('d0 = %.6f, d1 = %.6f, d2 = %.6f, d3 = %.6f\n', d0, d1, d2, d3);
fprintf('Mean Squared Error: E = %.6f\n\n', E_periodic);

% Plot the original data and the fitted linear + periodic model
figure;
plot(day, USDSEK, 'bo', 'DisplayName', 'Data');  % Original data points
hold on;
plot(day, g_periodic, 'r-', 'LineWidth', 2, 'DisplayName', 'Linear + Periodic Model');  % Fitted model
xlabel('Days');
ylabel('USD/SEK');
title('Linear + Periodic Model Fitting to Dollar Exchange Rate');
legend('Location', 'best');
grid on;
hold off;

%% Task c: Nonlinear Model Fitting using Gauss-Newton Method
% Task c: Improve the model by fitting the period L as a parameter using the Gauss-Newton method
% The model is now nonlinear in L: USDSEK = d0 + d1*day + d2*sin(2*pi*day/L) + d3*cos(2*pi*day/L)
% We will use the Gauss-Newton method to find the optimal parameters including L

fprintf('TASK C \n')

% Initialize parameters for nonlinear optimization
params = [d0, d1, d2, d3, L];  % Initial guesses for d0, d1, d2, d3, L
tol = 1e-10;  % Tolerance for convergence
max_iter = 100;  % Maximum number of iterations

% Begin Gauss-Newton iteration for nonlinear least squares fitting
for iter = 1:max_iter
    % Extract current parameter estimates
    d0 = params(1);
    d1 = params(2);
    d2 = params(3);
    d3 = params(4);
    L = params(5);
    
    % Compute the angle for the sine and cosine terms
    angle = 2 * pi * day / L;
    sin_angle = sin(angle);
    cos_angle = cos(angle);
    
    % Compute the model predictions using current parameters
    g_nonlinear = d0 + d1 * day + d2 * sin_angle + d3 * cos_angle;
    
    % Compute the residuals between observed data and model predictions
    residuals = USDSEK - g_nonlinear;
    
    % Compute the Jacobian matrix J with partial derivatives with respect to each parameter
    % The Jacobian matrix has dimensions N x 5, where N is the number of data points
    % Partial derivatives:
    %   - with respect to d0: 1
    %   - with respect to d1: day
    %   - with respect to d2: sin_angle
    %   - with respect to d3: cos_angle
    %   - with respect to L:
    %     ∂g/∂L = - (2πday / L^2) * (d2 * cos_angle - d3 * sin_angle)
    J = [ones(N, 1), ...                                % Partial derivative w.r.t d0
         day, ...                                       % Partial derivative w.r.t d1
         sin_angle, ...                                 % Partial derivative w.r.t d2
         cos_angle, ...                                 % Partial derivative w.r.t d3
         (-2 * pi * day .* (d2 * cos_angle - d3 * sin_angle)) / (L^2)];  % Partial derivative w.r.t L
         
    % Compute the Gauss-Newton update step delta
    delta = (J' * J) \ (J' * residuals);
    
    % Update the parameters
    params = params + delta';
    
    % Check for convergence
    if norm(delta) < tol
        fprintf('Gauss-Newton converged in %d iterations.\n\n', iter);
        break;
    end
end

% Check if the algorithm did not converge within the maximum iterations
if iter == max_iter
    warning('Gauss-Newton did not converge within the maximum number of iterations.');
end

% Extract the optimized parameters after convergence
d0 = params(1);
d1 = params(2);
d2 = params(3);
d3 = params(4);
L = params(5);

% Display the final optimized parameters
fprintf('TASK C - Nonlinear Model Coefficients:\n');
fprintf('d0 = %.6f\n', d0);
fprintf('d1 = %.6f\n', d1);
fprintf('d2 = %.6f\n', d2);
fprintf('d3 = %.6f\n', d3);
fprintf('L = %.6f\n', L);

% Compute the final fitted values using the optimized parameters
angle_final = 2 * pi * day / L; %Inside bracket in the formula g_final
g_final = d0 + d1 * day + d2 * sin(angle_final) + d3 * cos(angle_final);

% Compute the mean squared error (MSE) of the final model
E_final = mean((USDSEK - g_final).^2);
fprintf('Mean Squared Error: E = %.6f\n', E_final);

% Plot all models together for comparison
figure;
plot(day, USDSEK, 'bo', 'DisplayName', 'Data');  % Original data points
hold on;
plot(day, g_linear, 'r-', 'LineWidth', 2, 'DisplayName', 'Linear Model');  % Linear model
plot(day, g_periodic, 'k-', 'LineWidth', 2, 'DisplayName', 'Linear + Periodic Model');  % Linear + periodic model
plot(day, g_final, 'g-', 'LineWidth', 2, 'DisplayName', 'Nonlinear Model');  % Nonlinear model with optimized L
xlabel('Days');
ylabel('USD/SEK');
title('Comparison of Models Fitting to Dollar Exchange Rate');
legend('Location', 'best');
%grid on;
hold off;

