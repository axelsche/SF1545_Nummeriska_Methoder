clc;clear all; close all

% Parameters
L = 4.0;          % Length of the rod [m]
k = 2.2;          % Thermal conductivity [N/(K·s)]
t0 = 290.0;       % Temperature at x=0 [K]
t1 = 400.0;       % Temperature at x=L [K]
x_target = 3.0;   % Target position x = 3 m
h_initial = 0.2;  % Initial step size [m]
tolerance = 1e-5; % Desired accuracy (five correct decimal places)

% Initialize variables
h = h_initial;
T3_prev = 0;
error = inf;
iteration = 0;

% Arrays to store results for analysis
h_values = [];
T3_values = [];
errors = [];
orders = [];

fprintf('Iteration\tStep size h\tT(3)\t\tError\t\tEstimated Order p\n');

while error > tolerance
    iteration = iteration + 1;
    
    % Number of internal points
    n = round((L / h) - 1);
    
    % Generate x vector including boundary points
    x_full = linspace(0, L, n + 2)';  % Column vector from 0 to L with step size h
    xi = x_full(2:end-1);             % Internal points (excluding boundaries)
    
    % Compute Q(x) at internal points
    Q = 3000 * exp(-20 * (xi - 0.6 * L).^2) + 200;
    
    % Construct matrix A
    main_diag = 2 * ones(n, 1);           % Main diagonal elements
    off_diag = -1 * ones(n - 1, 1);       % Off-diagonal elements
    A = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);
    
    % Construct vector b
    b = (h^2 / k) * Q;
    
    % Adjust for boundary conditions
    b(1) = b(1) + t0;
    b(end) = b(end) + t1;
    
    % Solve the linear system A * T_internal = b
    T_internal = A \ b;
    
    % Include boundary temperatures
    T_full = [t0; T_internal; t1];
    
    % Find the approximation of T(3)
    index = find(abs(x_full - x_target) < 1e-10);  % Tolerance adjusted for floating-point errors
    
    if isempty(index)
        % Use interpolation if x_target is not in x_full
        T3 = interp1(x_full, T_full, x_target, 'linear');
    else
        T3 = T_full(index);
    end
    
    % Compute error if not the first iteration
    if iteration > 1
        error = abs(T3 - T3_prev);
        % Estimate order of accuracy
        p = log(abs(errors(end) / error)) / log(2);
        orders = [orders; p];
    else
        error = inf;
        p = NaN;
    end
    
    % Store results
    h_values = [h_values; h];
    T3_values = [T3_values; T3];
    errors = [errors; error];
    
    % Display results
    fprintf('%d\t\t%.6f\t%.6f\t%.6e\t%.2f\n', iteration, h, T3, error, p);
    
    % Prepare for next iteration
    T3_prev = T3;
    h = h / 2;  % Halve the step size
end

% Display final result
fprintf('\nFinal approximation of T(3) with five correct decimal places: %.10f K\n', T3);


% Plot the temperature distribution across the rod
figure;
plot(x_full, T_full, '-o', 'LineWidth', 1.5);
xline(x_target,'r--', LineWidth= 2)
xlabel('Position along the rod, x [m]');
ylabel('Temperature, T(x) [K]');
title('Temperature Distribution Along the Rod');
grid on;