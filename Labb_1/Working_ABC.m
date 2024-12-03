% Optimized MATLAB Code for Root Finding using Fixed-Point Iteration and Newton's Method

% ----------------------------
% Initialization
% ----------------------------
format long;          
clear variables;      
close ;               
clc;                  

% ----------------------------
% Define the Function and Its Derivative
% ----------------------------
f = @(x) x.^2 - 8*x - 12*sin(3*x + 1) + 19;       % Original function
df = @(x) 2*x - 8 - 36*cos(3*x + 1);            % Derivative of the function

% ----------------------------
% Plot the Function
% ----------------------------
x_plot = linspace(0, 7, 1000);                     % Define range for x
y_plot = f(x_plot);                                % Compute f(x) values

figure;
plot(x_plot, y_plot, 'b-', 'LineWidth', 1.5);      % Plot f(x) in blue solid line
hold on;
plot(x_plot, zeros(size(x_plot)), 'r--', 'LineWidth', 1.5); % Plot g(x)=0 in red dashed line
xlabel('x');
ylabel('f(x)');
title('Plot of the Function f(x) and g(x) = 0');
legend('f(x)', 'g(x) = 0');
grid on;
hold off;

% ----------------------------
% Ploting the derivative of the fixed point in a range of -1 < y < 1 to
% check for which points can be used for covergence
% ----------------------------


% Manually input the derivative of the fixed-point function
fixed_point_derivative = @(x) (2*x + 11 - 36*cos(3*x + 1)) / 19;

% Evaluate the derivative over the x range
derivative_values = fixed_point_derivative(x_plot);

% Plot the derivative
figure;
plot(x_plot, derivative_values, 'b', 'LineWidth', 2);
hold on;

% Plot horizontal lines at y = 1 and y = -1
y_limit = 1;
plot(x_plot, y_limit * ones(size(x_plot)), 'r--', 'LineWidth', 1);
plot(x_plot, -y_limit * ones(size(x_plot)), 'r--', 'LineWidth', 1);

% Add labels and title
xlabel('x');
ylabel('Derivative of fixed\_point(x)');
title('Derivative of Fixed-Point Function');

% Add legend
legend('Derivative', 'y = 1', 'y = -1', 'Location', 'Best');

% Enable grid
grid on;

% Set y-axis limits for better visualization (adjust as needed)
ylim([-5, 5]);

% Hold off plotting
hold off;


% ----------------------------
%% Fixed-Point Iteration Parameters
% ----------------------------
initial_guesses_fpi = [2, 2.5, 4, 4.5, 6, 6.5];    % Initial guesses for FPI
max_iter_fpi = 100;                               % Maximum iterations
tol_fpi = 1e-10;                                  % Tolerance for convergence

% Fixed-Point Iteration Formula
fixed_point = @(x) (x.^2 + 11*x - 12*sin(3*x + 1)) / 19 + 1;

% Preallocate Arrays
num_guesses = length(initial_guesses_fpi);
roots_fpi = zeros(1, num_guesses);                % To store roots
diff_fpi = cell(1, num_guesses);                  % To store differences

% ----------------------------
% Fixed-Point Iteration
% ----------------------------
fprintf('Fixed-Point Iteration Results:\n');
for j = 1:num_guesses
    x0 = initial_guesses_fpi(j);
    diffs = zeros(max_iter_fpi,1);                % Initialize differences
    for i = 1:max_iter_fpi
        x_new = fixed_point(x0);                  % Compute new x
        diffs(i) = abs(x_new - x0);               % Calculate difference
        
        % Check for Non-Convergence
        if i > 1 && diffs(i) / diffs(i-1) >= 1
            fprintf('Warning: Iteration %d for initial guess %.4f is not converging.\n', i, initial_guesses_fpi(j));
        end
        
        % Check for Convergence
        if diffs(i) < tol_fpi
            diffs = diffs(1:i);                    % Trim unused preallocated space
            break;
        end
        x0 = x_new;                                % Update x0 for next iteration
    end
    roots_fpi(j) = x0;                             % Store the found root
    diff_fpi{j} = diffs;                           % Store the differences
end

% Display Roots Found by Fixed-Point Iteration
fprintf('Roots found with Fixed-Point Iteration:\n');
disp(roots_fpi);

% ----------------------------
% Newton's Method Parameters
% ----------------------------
initial_guesses_newton = [2, 2.5, 4, 4.5, 6, 6.5]; % Initial guesses for Newton's Method
max_iter_newton = 100;                             % Maximum iterations
tol_newton = 1e-10;                                % Tolerance for convergence

% Preallocate Array
roots_newton = zeros(1, length(initial_guesses_newton));

% ----------------------------
% Newton's Method
% ----------------------------
fprintf('Newton''s Method Results:\n');
for j = 1:length(initial_guesses_newton)
    x = initial_guesses_newton(j);                   % Initial guess
    for k = 1:max_iter_newton
        fx = f(x);                                   % Compute f(x)
        dfx = df(x);                                 % Compute f'(x)
        
        % Check if derivative is too small
        if abs(dfx) < tol_newton
            warning('Derivative near zero at x = %.10f, may not converge.', x);
            break;
        end
        
        x_new = x - fx / dfx;                        % Update x using Newton's formula
        
        % Check for Convergence
        if abs(x_new - x) < tol_newton
            x = x_new;
            break;
        end
        x = x_new;                                    % Update x for next iteration
    end
    roots_newton(j) = x;                              % Store the found root
end

% Display Roots Found by Newton's Method
disp('Roots found with Newton''s Method:');
disp(roots_newton);
