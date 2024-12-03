
% ----------------------------
% Initialization
% ----------------------------
format long;              
clear variables;          
close all;                
clc;                      

% ----------------------------
% Define the Function and Its Derivative
% ----------------------------
f = @(x) x.^2 - 8*x - 12*sin(3*x + 1) + 19;        % Original function
df = @(x) 2*x - 8 - 36*cos(3*x + 1);             % Derivative of the function

% ----------------------------
% Reference Value (True Root)
% ----------------------------
x_ref = 1.972477260906544;    % Given reference root for error calculation, callculate using [Newton4Refferensing.m]

% ----------------------------
% Initial Parameters
% ----------------------------
x0 = 2;                        % Initial guess
max_iter = 100;                % Maximum number of iterations
tol = 1e-15;                   % Tolerance for convergence

% ----------------------------
% Initialize Error Arrays
% ----------------------------
errors_fpi = zeros(max_iter, 1);        % Errors for Fixed-Point Iteration
errors_newton = zeros(max_iter, 1);     % Errors for Newton's Method

% ----------------------------
% Fixed-Point Iteration
% ----------------------------
% Fixed-Point Iteration Formula: x_new = g(x) = (x^2 + 11x - 12*sin(3x + 1))/19 + 1
fixed_point = @(x) (x.^2 + 11*x - 12*sin(3*x + 1)) / 19 + 1;

x_fpi = x0;                     % Initialize FPI with initial guess
iter_fpi = 0;                   % Iteration counter for FPI

for i = 1:max_iter
    x_new = fixed_point(x_fpi);          % Compute new x using FPI formula
    errors_fpi(i) = abs(x_new - x_ref);  % Calculate error w.r.t reference root
    
    % Check for convergence
    if abs(x_new - x_fpi) < tol
        iter_fpi = i;                     % Record number of iterations
        break;
    end
    
    % Update for next iteration
    x_fpi = x_new;
end

% Trim the unused portion of the error array
if iter_fpi == 0
    iter_fpi = max_iter;
end
errors_fpi = errors_fpi(1:iter_fpi);

% ----------------------------
% Newton's Method
% ----------------------------
x_newton = x0;                  % Initialize Newton's Method with initial guess
iter_newton = 0;                % Iteration counter for Newton's Method

for i = 1:max_iter
    fx = f(x_newton);           % Compute f(x)
    dfx = df(x_newton);         % Compute f'(x)
    
    % Check if derivative is too small
    if abs(dfx) < tol
        warning('Derivative near zero at iteration %d, x = %.10f. Method may not converge.', i, x_newton);
        break;
    end
    
    x_new = x_newton - fx / dfx;           % Update x using Newton's formula
    errors_newton(i) = abs(x_new - x_ref); % Calculate error w.r.t reference root
    
    % Check for convergence
    if abs(x_new - x_newton) < tol
        iter_newton = i;                   % Record number of iterations
        break;
    end
    
    % Update for next iteration
    x_newton = x_new;
end

% Trim the unused portion of the error array
if iter_newton == 0
    iter_newton = max_iter;
end
errors_newton = errors_newton(1:iter_newton);

% ----------------------------
% Plotting the Errors vs. Number of Iterations
% ----------------------------
figure;
semilogy(1:iter_fpi, errors_fpi, 'r-', 'LineWidth', 2); hold on;    % FPI Error Line
semilogy(1:iter_fpi, errors_fpi, 'ro', 'MarkerFaceColor', 'r');    % FPI Error Points
semilogy(1:iter_newton, errors_newton, 'b--', 'LineWidth', 2);     % Newton's Method Error Line
semilogy(1:iter_newton, errors_newton, 'bs', 'MarkerFaceColor', 'b');% Newton's Method Error Points
hold off;

% Adding labels and title
xlabel('Number of Iterations');
ylabel('Error |x_n - x^*|');
title('Error vs. Number of Iterations for FPI and Newton''s Method');
legend('Fixed-Point Iteration', 'FPI Points', 'Newton''s Method', 'Newton''s Points', 'Location', 'best');
grid on;

% ----------------------------
% Log-Log Plot for Convergence Order
% ----------------------------
figure;
loglog(errors_fpi(1:end-1), errors_fpi(2:end), 'r-', 'LineWidth', 2); hold on; % FPI Convergence Line
loglog(errors_fpi(1:end-1), errors_fpi(2:end), 'ro', 'MarkerFaceColor', 'r');     % FPI Convergence Points
loglog(errors_newton(1:end-1), errors_newton(2:end), 'b--', 'LineWidth', 2);      % Newton's Method Convergence Line
loglog(errors_newton(1:end-1), errors_newton(2:end), 'bs', 'MarkerFaceColor', 'b');% Newton's Method Convergence Points
hold off;

% labels and title
xlabel('Error |e_n|');
ylabel('Error |e_{n+1}|');
title('Log-Log Plot of Error for FPI and Newton''s Method');
legend('Fixed-Point Iteration', 'FPI Points', 'Newton''s Method', 'Newton''s Points', 'Location', 'best');
grid on;

% ----------------------------
% Display Results
% ----------------------------
fprintf('Fixed-Point Iteration:\n');
fprintf('Root found: %.12f\n', x_fpi);
fprintf('Number of iterations: %d\n', iter_fpi);
fprintf('Final error: %.12e\n\n', errors_fpi(end));

fprintf('Newton''s Method:\n');
fprintf('Root found: %.12f\n', x_newton);
fprintf('Number of iterations: %d\n', iter_newton);
fprintf('Final error: %.12e\n', errors_newton(end));
