%% 2b Print two circles and the string between them
A = [-1.5, 3.0];
B = [1.0, 1.0];

xa = A(1);
ya = A(2);

xb = B(1);
yb = B(2);

ra = 1.5;
rb = 0.8;

% Initial guesses for Newton's method
x1 = -5; y1 = 10;
x2 = 0; y2 = 10;

[n_vector, niter] = Newton(@F, @J, x1, y1, x2, y2, xa, xb, ya, yb, ra, rb);
fprintf('Solution:\n');
fprintf('(%.10f, %.10f), (%.10f, %.10f)\n', n_vector(1), n_vector(2), n_vector(3), n_vector(4));
fprintf('Number of iterations: %d\n', niter)

figure(1)
hold on;
plot_circle(xa, ya, ra, 'b');
plot_circle(xb, yb, rb, 'r');

plot([n_vector(1) n_vector(3)], [n_vector(2) n_vector(4)], 'k', 'LineWidth', 1.5);
plot(n_vector(1), n_vector(2), 'k.');
plot(n_vector(3), n_vector(4), 'k.');
axis equal;
xlabel('X-axis');
ylabel('Y-axis');
hold off;

%% 2c Calculate the length of the string around the three circles

% Coordinates of the centers
A = [-1, 1.5];
B = [3.0, 0.5];
C = [0.0, -2.0];

% Circle radii
ra = 1;
rb = 1.2;
rc = 1.7;

% Initial guesses for the six points
x1 = 0;   y1 = 3;
x2 = 4;   y2 = 2;
x3 = -3;  y3 = 2;
x4 = -3;  y4 = -2;
x5 = 2;   y5 = -3;
x6 = 4;   y6 = 0;

% Use Newton's method to find tangent points
[newton_vector1, niter1] = Newton(@F, @J, x1, y1, x2, y2, A(1), B(1), A(2), B(2), ra, rb);
[newton_vector2, niter2] = Newton(@F, @J, x3, y3, x4, y4, A(1), C(1), A(2), C(2), ra, rc);
[newton_vector3, niter3] = Newton(@F, @J, x5, y5, x6, y6, C(1), B(1), C(2), B(2), rc, rb);

% Update the six points with results from Newton's method
x1 = newton_vector1(1); y1 = newton_vector1(2);
x2 = newton_vector1(3); y2 = newton_vector1(4);
x3 = newton_vector2(1); y3 = newton_vector2(2);
x4 = newton_vector2(3); y4 = newton_vector2(4);
x5 = newton_vector3(1); y5 = newton_vector3(2);
x6 = newton_vector3(3); y6 = newton_vector3(4);

% Create a figure
figure(2)
hold on;

% Plot circles
plot_circle(A(1), A(2), ra, 'b');
plot_circle(B(1), B(2), rb, 'g');
plot_circle(C(1), C(2), rc, 'r');

% Plot lines
plot([x1 x2], [y1 y2], 'k', 'LineWidth', 1.5);
plot([x3 x4], [y3 y4], 'k', 'LineWidth', 1.5);
plot([x5 x6], [y5 y6], 'k', 'LineWidth', 1.5);

% Plot points
plot(x1, y1, 'bo');
plot(x2, y2, 'go');
plot(x3, y3, 'bo');
plot(x4, y4, 'ro');
plot(x5, y5, 'ro');
plot(x6, y6, 'go');

axis equal;
xlabel('X-axis');
ylabel('Y-axis');
hold off;

% Call function to calculate the length of the string
length = string_length(x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, ra, rb, rc);
fprintf('Length of the string: %.8f \n', length)

%% 2d Calculate the uncertainty in the string's length

f = 0.01;  % Define the uncertainty in the coordinates

% Original centers
A = [-1, 1.5];
B = [3.0, 0.5];
C = [0.0, -2.0];

% Original string length without uncertainty
true_length = lengthO(A, B, C, ra, rb, rc);

% Create variations of the centers with uncertainty
Af1 = A + [f, 0];  % Move A slightly in the x-direction
Af2 = A + [0, f];  % Move A slightly in the y-direction
Bf1 = B + [f, 0];  % Move B slightly in the x-direction
Bf2 = B + [0, f];  % Move B slightly in the y-direction
Cf1 = C + [f, 0];  % Move C slightly in the x-direction
Cf2 = C + [0, f];  % Move C slightly in the y-direction
raf = ra + f;
rbf = rb + f;
rcf = rc + f;

% Calculate the string length for various uncertain versions of the centers and radii
uncertain_lengths = [
    lengthO(Af1, B, C, ra, rb, rc);
    lengthO(Af2, B, C, ra, rb, rc);
    lengthO(A, Bf1, C, ra, rb, rc);
    lengthO(A, Bf2, C, ra, rb, rc);
    lengthO(A, B, Cf1, ra, rb, rc);
    lengthO(A, B, Cf2, ra, rb, rc);
    lengthO(A, B, C, raf, rb, rc);
    lengthO(A, B, C, ra, rbf, rc);
    lengthO(A, B, C, ra, rb, rcf);
];

% Calculate the differences between the uncertain lengths and the true length
diff_lengths = abs(uncertain_lengths - true_length);

% Calculate total uncertainty as the sum of differences
total_uncertainty = sum(diff_lengths);

% Alternatively, take the maximum difference as uncertainty
max_uncertainty = max(diff_lengths);

% Print uncertainties
fprintf('Total uncertainty (sum of differences): %.7f \n', total_uncertainty)
fprintf('Maximum uncertainty: %.7f \n', max_uncertainty)

%% Functions

function plot_circle(x_m, y_m, r, color)
    theta = linspace(0, 2*pi, 100);
    x = x_m + r * cos(theta);
    y = y_m + r * sin(theta);
    plot(x, y, color, 'LineWidth', 2)
end

% Function to calculate the straight-line distance between two points
function l = l_straight(x1, y1, x2, y2)
    l = sqrt((x2 - x1)^2 + (y2 - y1)^2);
end

% Function to calculate the arc length between two points on a circle
function l = l_circle(x1, y1, x2, y2, r)
    chord_length = l_straight(x1, y1, x2, y2);
    theta = 2 * asin(chord_length / (2 * r));
    l = theta * r;
end

% Function that calculates the length of the string
function length = string_length(x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, ra, rb, rc)
    length = l_straight(x1, y1, x2, y2) ...
           + l_straight(x3, y3, x4, y4) ...
           + l_straight(x5, y5, x6, y6) ...
           + l_circle(x1, y1, x3, y3, ra) ...
           + l_circle(x2, y2, x6, y6, rb) ...
           + l_circle(x5, y5, x4, y4, rc);
end

% Newton's method to solve the system of equations
function [n_vector, niter] = Newton(F, J, x1, y1, x2, y2, xa, xb, ya, yb, ra, rb)
    tol = 1e-10;      % Tolerance for convergence (stopping criterion)
    max_iter = 100;   % Max number of iterations for Newton's method
    initial_vector = [x1; y1; x2; y2];
    niter = 0;
    res = tol + 1;
    while (res > tol) && (niter < max_iter)
        J_Newton = J(initial_vector(1), initial_vector(2), initial_vector(3), initial_vector(4), xa, xb, ya, yb);
        F_Newton = F(initial_vector(1), initial_vector(2), initial_vector(3), initial_vector(4), xa, xb, ya, yb, ra, rb);
        delta = J_Newton \ F_Newton;
        n_vector = initial_vector - delta;
        res = norm(delta);
        initial_vector = n_vector;
        niter = niter + 1;
    end
end

% Define the function F
function result = F(x1, y1, x2, y2, xa, xb, ya, yb, ra, rb)
    result = [
        (x1 - xa)^2 + (y1 - ya)^2 - ra^2;
        (x2 - xb)^2 + (y2 - yb)^2 - rb^2;
        (x1 - x2) * (x1 - xa) + (y1 - y2) * (y1 - ya);
        (x1 - x2) * (x2 - xb) + (y1 - y2) * (y2 - yb)
    ];
end

% Define the Jacobian matrix J
function result = J(x1, y1, x2, y2, xa, xb, ya, yb)
    result = [
        2 * (x1 - xa), 2 * (y1 - ya), 0, 0;
        0, 0, 2 * (x2 - xb), 2 * (y2 - yb);
        (2 * x1 - x2 - xa), (2 * y1 - y2 - ya), (-x1 + xa), (-y1 + ya);
        (-x2 + xb), (-y2 + yb), (-2 * x2 + x1 + xb), (-2 * y2 + y1 + yb)
    ];
end

% Function to calculate the length of the string with variations in the centers
function uncertain_length = lengthO(A, B, C, ra, rb, rc)
    % Initial guesses for the six points
    x1 = 0;   y1 = 3;
    x2 = 4;   y2 = 2;
    x3 = -3;  y3 = 2;
    x4 = -3;  y4 = -2;
    x5 = 2;   y5 = -3;
    x6 = 4;   y6 = 0;

    % Use Newton's method to find exact coordinates with the uncertain centers
    [n_vector1, ~] = Newton(@F, @J, x1, y1, x2, y2, A(1), B(1), A(2), B(2), ra, rb);
    [n_vector2, ~] = Newton(@F, @J, x3, y3, x4, y4, A(1), C(1), A(2), C(2), ra, rc);
    [n_vector3, ~] = Newton(@F, @J, x5, y5, x6, y6, C(1), B(1), C(2), B(2), rc, rb);

    % Update the six points with results from Newton's method
    x1 = n_vector1(1); y1 = n_vector1(2);
    x2 = n_vector1(3); y2 = n_vector1(4);
    x3 = n_vector2(1); y3 = n_vector2(2);
    x4 = n_vector2(3); y4 = n_vector2(4);
    x5 = n_vector3(1); y5 = n_vector3(2);
    x6 = n_vector3(3); y6 = n_vector3(4);

    % Calculate the string's length with the updated points
    uncertain_length = string_length(x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, ra, rb, rc);
end

