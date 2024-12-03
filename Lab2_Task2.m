 % Clear the workspace and close all figures
clc; clear ; close ;

%% Lab 2.a
disp('-------------------------Part_A-------------------------');

% Parameters for the damped pendulum
m = 0.6;          % Mass (kg)
mu = 0.2;         % Damping coefficient
g = 9.81;         % Acceleration due to gravity (m/s^2)
L = 1.5;          % Length of the pendulum (m)
h = 0.025;        % Time step size
t_end = 5;        % End time (s)
phi0 = 0.5;       % Initial angle (radians)
phi_dot0 = 0;     % Initial angular velocity (rad/s)

% Define the ODE system:
% y(1) = φ (angle)
% y(2) = φ' (angular velocity)
% y(1)' = y(2)
% y(2)' = - (μ/m) * y(2) - (g/L) * sin(y(1))
f = @(y) [y(2); - (mu/m) * y(2) - (g/L) * sin(y(1))];

% Time vector
t = 0:h:t_end;
N = length(t);

% Initialize solution arrays
y_euler = zeros(2, N);  % For Forward Euler method
y_rk4 = zeros(2, N);    % For Runge-Kutta 4 method

% Set initial conditions
y_euler(:, 1) = [phi0; phi_dot0];
y_rk4(:, 1) = [phi0; phi_dot0];

% Forward Euler method
for n = 1:N-1
    y_euler(:, n+1) = y_euler(:, n) + h * f(y_euler(:, n));
end

% Runge-Kutta 4 method
for n = 1:N-1
    % Compute the Runge-Kutta 4 coefficients
    k1 = h * f(y_rk4(:, n));
    k2 = h * f(y_rk4(:, n) + 0.5 * k1);
    k3 = h * f(y_rk4(:, n) + 0.5 * k2);
    k4 = h * f(y_rk4(:, n) + k3);
    
    % Update the solution
    y_rk4(:, n+1) = y_rk4(:, n) + (k1 + 2*k2 + 2*k3 + k4)/6;
end

% Plot the solutions
figure;
plot(t, y_euler(1, :), 'b-', 'DisplayName', 'Forward Euler');
hold on;
plot(t, y_rk4(1, :), 'r-', 'DisplayName', 'Runge-Kutta 4');
xlabel('Time t (s)');
ylabel('\phi(t) (rad)');
title('Solution for damped pendulum using Forward Euler and Runge-Kutta 4 methods');
legend;
hold off;

% Display the approximations of φ(5) for both methods
fprintf('Forward Euler: φ(5) ≈ %.5f\n', y_euler(1, end));
fprintf('Runge-Kutta 4: φ(5) ≈ %.5f\n', y_rk4(1, end));


%% Lab 2.b
disp('-------------------------part_B-------------------------');

% Error tolerances
tol_euler = 0.5e-2;   % Tolerance for Forward Euler method
tol_rk4 = 0.5e-6;     % Tolerance for Runge-Kutta 4 method

% Arrays to store approximations and step sizes for error analysis
phi_euler_results = [];
h_values_euler = [];

phi_rk4_results = [];
h_values_rk4 = [];

% Initial step sizes
h_euler = 0.1;    % Starting step size for Euler method
h_rk4 = 0.05;     % Starting step size for RK4 method

% Initial error values (set to large values to enter the loops)
error_euler = inf;
error_rk4 = inf;

% Previous approximations (initialize arbitrarily)
prev_phi_euler = 0;
prev_phi_rk4 = 0;

% Calculations with Forward Euler method, refining h until error is below tolerance
disp(' ')
fprintf('Calculations with Forward Euler method:\n');
while error_euler > tol_euler
    [t_euler, y_euler] = forward_euler(f, h_euler, [phi0; phi_dot0], t_end);
    current_phi_euler = y_euler(1, end);
    error_euler = abs(current_phi_euler - prev_phi_euler);
    
    % Store results for error analysis
    phi_euler_results = [phi_euler_results, current_phi_euler];
    h_values_euler = [h_values_euler, h_euler];
    
    fprintf('h = %.5f, φ(5) ≈ %.10f, Error = %.10f\n', h_euler, current_phi_euler, error_euler);
    
    % Update for next iteration
    prev_phi_euler = current_phi_euler;
    h_euler = h_euler / 2;
end

% Calculations with Runge-Kutta 4 method, refining h until error is below tolerance
disp(' ')
fprintf('Calculations with Runge-Kutta 4 method:\n');
prev_phi_rk4 = 0;

while error_rk4 > tol_rk4
    [t_rk4, y_rk4] = runge_kutta_4(f, h_rk4, [phi0; phi_dot0], t_end);
    current_phi_rk4 = y_rk4(1, end);
    error_rk4 = abs(current_phi_rk4 - prev_phi_rk4);
    
    % Store results for error analysis
    phi_rk4_results = [phi_rk4_results, current_phi_rk4];
    h_values_rk4 = [h_values_rk4, h_rk4];
    
    fprintf('h = %.5f, φ(5) ≈ %.10f, Error = %.10f\n', h_rk4, current_phi_rk4, error_rk4);
    
    % Update for next iteration
    prev_phi_rk4 = current_phi_rk4;
    h_rk4 = h_rk4 / 2;
end

% Estimate order of accuracy for Forward Euler method
disp(' ')
fprintf('Order of accuracy for Forward Euler method:\n');
for i = 1:length(h_values_euler) - 2
    e1 = phi_euler_results(i) - phi_euler_results(i+1);
    e2 = phi_euler_results(i+1) - phi_euler_results(i+2);
    ratio = e1 / e2;
    order = log2(abs(ratio));
    fprintf('h = %.5f, Ratio = %.5f, Estimated order = %.2f\n', h_values_euler(i+2), ratio, order);
end

% Estimate order of accuracy for Runge-Kutta 4 method
disp(' ')
fprintf('Order of accuracy for Runge-Kutta 4 method:\n');
for i = 1:length(h_values_rk4) - 2
    e1 = phi_rk4_results(i) - phi_rk4_results(i+1);
    e2 = phi_rk4_results(i+1) - phi_rk4_results(i+2);
    ratio = e1 / e2;
    order = log2(abs(ratio));
    fprintf('h = %.5f, Ratio = %.5f, Estimated order = %.2f\n', h_values_rk4(i+1), ratio, order);
end

% Display the final approximations
disp(' ')
fprintf('Results:\n');
fprintf('Forward Euler approximation for φ(5) with two correct decimals: φ(5) ≈ %.10f\n', current_phi_euler);
fprintf('Runge-Kutta 4 approximation for φ(5) with six correct decimals: φ(5) ≈ %.10f\n', current_phi_rk4);

% --- Function Definitions ---

% Function for the Forward Euler method
function [t, y] = forward_euler(f, h, y0, t_end)
    % Forward Euler method to solve y' = f(y)
    %
    % Inputs:
    %   f     - Function handle representing the ODE system
    %   h     - Step size
    %   y0    - Initial condition vector
    %   t_end - End time
    %
    % Outputs:
    %   t - Time vector
    %   y - Solution matrix, each column corresponds to time step

    t = 0:h:t_end;
    N = length(t);
    y = zeros(length(y0), N);
    y(:, 1) = y0;
    for n = 1:N-1
        y(:, n+1) = y(:, n) + h * f(y(:, n));
    end
end

% Function for the Runge-Kutta 4 method
function [t, y] = runge_kutta_4(f, h, y0, t_end)
    % Runge-Kutta 4 method to solve y' = f(y)
    %
    % Inputs:
    %   f     - Function handle representing the ODE system
    %   h     - Step size
    %   y0    - Initial condition vector
    %   t_end - End time
    %
    % Outputs:
    %   t - Time vector
    %   y - Solution matrix, each column corresponds to time step

    t = 0:h:t_end;
    N = length(t);
    y = zeros(length(y0), N);
    y(:, 1) = y0;
    for n = 1:N-1
        % Compute the Runge-Kutta 4 coefficients
        k1 = h * f(y(:, n));
        k2 = h * f(y(:, n) + 0.5 * k1);
        k3 = h * f(y(:, n) + 0.5 * k2);
        k4 = h * f(y(:, n) + k3);
        
        % Update the solution
        y(:, n+1) = y(:, n) + (k1 + 2*k2 + 2*k3 + k4)/6;
    end
end
