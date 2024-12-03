% ----------------------------
% Define the Function and Its Derivative
% ----------------------------
f = @(x) x.^2 - 8*x - 12*sin(3*x + 1) + 19;       % Original function
df = @(x) 2*x - 8 - 36*cos(3*x + 1);              % Derivative of the function

% ----------------------------
% Newton's Method Parameters for Reference Solution
% ----------------------------
initial_guess = 2;                                % Initial guess for Newton's method
max_iter_newton_ref = 1000;                       % Maximum iterations
accuracy = 15;                                    % Number of correct decimal digits
tol_newton_ref = 10^(-accuracy);                  % Tolerance based on desired accuracy

% ----------------------------
% Compute a very accurate root x* using Newton's Method
% ----------------------------
x_star = initial_guess;                           % Start with the initial guess
for k = 1:max_iter_newton_ref
    fx = f(x_star);                               % Compute f(x)
    dfx = df(x_star);                             % Compute f'(x)
    
    % Check if the derivative is close to zero (risk of non-convergence)
    if abs(dfx) < tol_newton_ref
        warning('Derivative near zero at x = %.10f, may not converge.', x_star);
        break;
    end
    
    % Newton's update step
    x_new = x_star - fx / dfx;
    
    % Check for convergence
    if abs(x_new - x_star) < tol_newton_ref
        x_star = x_new;
        break;
    end
    
    % Update x* for the next iteration
    x_star = x_new;
end

% Display the highly accurate reference solution
fprintf('Accurate reference solution (x*) with %d decimal digits: %.15f\n', accuracy, x_star);
