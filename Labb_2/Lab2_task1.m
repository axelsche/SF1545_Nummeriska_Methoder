% Main script for parts (a), (b), and (c)
clear ; close ; clc

% Parameters
R = 3; 
L = 3 * sqrt(2);
tolerance_a = 0.5 * 1e-4; % Error tolerance for part (a)
tolerance_b = 1e-3; % Error tolerance for part (b)
initial_n = 30; % Starting with n = 30 for both parts

% Define g(r) function
g = @(r) (3 * r.^3 .* exp(-r)) ./ (1 + (1/3) * sin(8 * r / 5));

% -------------------------------------------------
% Part (a): Compute volume using polar coordinates
% -------------------------------------------------
disp('---------------------------------Part_A---------------------------------')
disp('')

V0 = g(R) * R^2 * pi;

% Trapezoidal Rule for part (a)
n_a = initial_n;
V_approx_trap_a = 0;
error_trap_a = inf;
while error_trap_a > tolerance_a
    integral_trap_a = trapezoidal_1D(@(r) g(r) .* r, 0, R, n_a);
    V_new_trap_a = V0 - 2 * pi * integral_trap_a;
    error_trap_a = abs(V_new_trap_a - V_approx_trap_a);
    fprintf('Part (a) Trapezoidal: n = %d, V = %.8f, Error = %.8f\n', n_a, V_new_trap_a, error_trap_a);
    
    V_approx_trap_a = V_new_trap_a;
    n_a = n_a * 2; % Double the number of intervals
end

disp(' ')

% Simpson's Rule for part (a)
n_a = initial_n;
V_approx_simpson_a = 0;
error_simpson_a = inf;
while error_simpson_a > tolerance_a
    % Ensure n_a is even for Simpson's rule
    if mod(n_a, 2) ~= 0
        n_a = n_a + 1;
    end
    integral_simpson_a = simpson_1D(@(r) g(r) .* r, 0, R, n_a);
    V_new_simpson_a = V0 - 2 * pi * integral_simpson_a;
    error_simpson_a = abs(V_new_simpson_a - V_approx_simpson_a);
    fprintf('Part (a) Simpson: n = %d, V = %.8f, Error = %.8f\n', n_a, V_new_simpson_a, error_simpson_a);
    
    V_approx_simpson_a = V_new_simpson_a;
    n_a = n_a * 2; % Double the number of intervals
end

disp(' ')
fprintf('\n========== FINAL VOLUME (TRAPEZOIDAL RULE, PART A): %.8f ==========\n', V_approx_trap_a);
fprintf('========== FINAL VOLUME (SIMPSONS RULE, PART A): %.8f ==========\n\n', V_approx_simpson_a);


% -------------------------------------------------
% Part (b): Compute volume using square opening
% -------------------------------------------------
disp('---------------------------------Part_B---------------------------------')
disp('')

% Calculate g(R) as a constant for part (b)
g_R = g(R);

% Initialize variables for iteration in part (b)
n_b = initial_n;
V_approx_b = 0;
error_b = inf;

% Iteratively refine n until the error is within tolerance
while error_b > tolerance_b
    % Calculate the volume using the 2D trapezoidal rule
    V_new_b = trapezoidal_2D_square(@(x, y) g_R - g(sqrt(x.^2 + y.^2)), L, n_b);
    
    % Calculate the error as the difference between successive approximations
    error_b = abs(V_new_b - V_approx_b);
    fprintf('Part (b) 2D Trapezoidal: n = %d, V = %.8f, Error = %.8f\n', n_b, V_new_b, error_b);
    
    % Update the approximation and double n for the next iteration
    V_approx_b = V_new_b;
    n_b = n_b * 2;
end

disp(' ')
fprintf('\n\n===== FINAL VOLUME APPROXIMATION (PART B): %.8f =====\n\n', V_approx_b);

% -------------------------------------------------
% Part (c): Error Analysis and Convergence Rate
% -------------------------------------------------
disp('---------------------------------Part_C---------------------------------')
disp('')
% Define a large n for an accurate reference solution
reference_n = 2048;
% Ensure reference_n is even for Simpson's rule
if mod(reference_n, 2) ~= 0
    reference_n = reference_n + 1;
end
V_ref_trap = V0 - 2 * pi * trapezoidal_1D(@(r) g(r) .* r, 0, R, reference_n);
V_ref_simpson = V0 - 2 * pi * simpson_1D(@(r) g(r) .* r, 0, R, reference_n);
V_ref_trap_2D = trapezoidal_2D_square(@(x, y) g_R - g(sqrt(x.^2 + y.^2)), L, reference_n);

% Define h values for analysis
h_values = L ./ [30, 60, 120, 240];

% Initialize error and convergence order arrays
errors_trap = zeros(size(h_values));
errors_simpson = zeros(size(h_values));
errors_trap_2D = zeros(size(h_values));
orders_trap = zeros(size(h_values));
orders_simpson = zeros(size(h_values));
orders_trap_2D = zeros(size(h_values));

% Loop over h values and calculate errors for each method
for i = 1:length(h_values)
    h = h_values(i);
    n = round(L / h);
    
    % Ensure n is even for Simpson's rule
    if mod(n, 2) ~= 0
        n = n + 1;
    end
    
    % 1D Trapezoidal Rule Error
    V_trap = V0 - 2 * pi * trapezoidal_1D(@(r) g(r) .* r, 0, R, n);
    E_trap = abs(V_ref_trap - V_trap);
    errors_trap(i) = E_trap;
    
    % 1D Simpson's Rule Error
    V_simpson = V0 - 2 * pi * simpson_1D(@(r) g(r) .* r, 0, R, n);
    E_simpson = abs(V_ref_simpson - V_simpson);
    errors_simpson(i) = E_simpson;
    
    % 2D Trapezoidal Rule Error
    V_trap_2D = trapezoidal_2D_square(@(x, y) g_R - g(sqrt(x.^2 + y.^2)), L, n);
    E_trap_2D = abs(V_ref_trap_2D - V_trap_2D);
    errors_trap_2D(i) = E_trap_2D;
    
    % Compute convergence orders after the first iteration
    if i > 1
        h_prev = h_values(i-1);
        E_trap_prev = errors_trap(i-1);
        E_simpson_prev = errors_simpson(i-1);
        E_trap_2D_prev = errors_trap_2D(i-1);
        
        order_trap = log(E_trap / E_trap_prev) / log(h / h_prev);
        orders_trap(i) = order_trap;
        
        order_simpson = log(E_simpson / E_simpson_prev) / log(h / h_prev);
        orders_simpson(i) = order_simpson;
        
        order_trap_2D = log(E_trap_2D / E_trap_2D_prev) / log(h / h_prev);
        orders_trap_2D(i) = order_trap_2D;
    end
end

% Display convergence orders
fprintf('\nConvergence Orders:\n');
for i = 2:length(h_values)
    fprintf('h = %.5f: Trapezoidal Order = %.2f, Simpson Order = %.2f, 2D Trapezoidal Order = %.2f\n', ...
        h_values(i), orders_trap(i), orders_simpson(i), orders_trap_2D(i));
end

% Plot error vs h on a log-log scale
figure;
loglog(h_values, errors_trap, '-o', 'DisplayName', '1D Trapezoidal');
hold on;
loglog(h_values, errors_simpson, '-s', 'DisplayName', '1D Simpson');
loglog(h_values, errors_trap_2D, '-^', 'DisplayName', '2D Trapezoidal');

% Add reference lines for expected convergence rates
h_ref = h_values;
C_trap = errors_trap(1) / h_values(1)^2;
C_simpson = errors_simpson(1) / h_values(1)^4;
C_trap_2D = errors_trap_2D(1) / h_values(1)^2;
loglog(h_ref, C_trap * h_ref.^2, '--', 'DisplayName', 'O(h^2)');
loglog(h_ref, C_simpson * h_ref.^4, '--', 'DisplayName', 'O(h^4)');
loglog(h_ref, C_trap_2D * h_ref.^2, '--', 'DisplayName', 'O(h^2)');

xlabel('Step size h');
ylabel('Absolute Error');
title('Error vs. Step Size for Different Methods');
legend('Location', 'best');
grid on;

% Functions used in the script
% -----------------------------

% Trapezoidal Rule for 1D integration
function I = trapezoidal_1D(f, a, b, n)
    h = (b - a) / n;
    x = linspace(a, b, n+1);
    fx = f(x);
    I = (h/2) * (fx(1) + 2*sum(fx(2:end-1)) + fx(end));
end

% Simpson's Rule for 1D integration
function I = simpson_1D(f, a, b, n)
    h = (b - a) / n;
    x = linspace(a, b, n+1);
    fx = f(x);
    I = (h/3) * (fx(1) + 4*sum(fx(2:2:end-1)) + 2*sum(fx(3:2:end-2)) + fx(end));
end

% Trapezoidal Rule for 2D integration over a square domain
function I = trapezoidal_2D_square(f, L, n)
    h = L / n;
    x = linspace(-L/2, L/2, n+1);
    y = linspace(-L/2, L/2, n+1);
    [X, Y] = meshgrid(x, y);
    fx = f(X, Y);
    % Apply trapezoidal weights
    fx([1 end], :) = fx([1 end], :) / 2;
    fx(:, [1 end]) = fx(:, [1 end]) / 2;
    I = h^2 * sum(fx(:));
end
