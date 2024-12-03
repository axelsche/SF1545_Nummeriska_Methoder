% Rensa arbetsytan och figurer
clear all;
close all;
clc;

% Parametrar
m = 0.75;               % Massa [kg]
g = 9.81;               % Tyngdacceleration [m/s^2]
Kxx = 0.004;            % Luftmotståndskonstant i x-led [kg/m]
Kyy = 0.078;            % Luftmotståndskonstant i y-led [kg/m]
k1 = Kxx / m;           % Luftmotstånd per massenhet i x-led [1/m]
k2 = Kyy / m;           % Luftmotstånd per massenhet i y-led [1/m]
v0 = 18.5;              % Starthastighet [m/s]
alpha = pi / 4;         % Kastvinkel [rad]
h0 = 0.5;               % Start höjd [m]

% Begynnelsevärden
vx0 = v0 * cos(alpha);  % Initial hastighet i x-led [m/s]
vy0 = v0 * sin(alpha);  % Initial hastighet i y-led [m/s]
x0 = 0;                 % Initial x-position [m]
y0 = h0;                % Initial y-position [m]

% Tolerans för felet
tolerance = 1e-6;

% Initial steglängd
h_initial = 0.01;

% Beräkna banan och nedslagsplatsen
[x_hit, error_estimate, X1] = compute_trajectory(m, g, k1, k2, x0, y0, vx0, vy0, h_initial, tolerance);

% Visa resultatet
fprintf('Nedslagsplatsens x-koordinat: %.8f meter\n', x_hit);
fprintf('Feluppskattning: %.2e\n', error_estimate);

% Visualisering av banan
figure;
plot(X1(:,1), X1(:,2), 'b-', 'LineWidth', 1.5);
hold on;
plot(x_hit, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');  % Markera nedslagsplatsen
xlabel('x [m]');
ylabel('y [m]');
title('Varpastenens bana');
grid on;
legend('Bana', 'Nedslagsplats', 'Location', 'Best');

function [x_hit, error_estimate, X1] = compute_trajectory(m, g, k1, k2, x0, y0, vx0, vy0, h_initial, tolerance)
    % Initial steglängd
    h_current = h_initial;
    % Initialisera feluppskattning
    error_estimate = Inf;
    
    while error_estimate > tolerance
        % Lös med steglängd h_current
        [t1, X1] = RK4_system(h_current, x0, y0, vx0, vy0, k1, k2, g);
        x1_hit = interpolate_landing(X1);
        
        % Lös med dubbel steglängd
        [t2, X2] = RK4_system(2*h_current, x0, y0, vx0, vy0, k1, k2, g);
        x2_hit = interpolate_landing(X2);
        
        % Beräkna feluppskattning
        error_estimate = abs(x1_hit - x2_hit);
        
        if error_estimate > tolerance
            % Halvera steglängden
            h_current = h_current / 2;
        else
            % Tillräcklig noggrannhet uppnådd
            break;
        end
    end
    
    x_hit = x1_hit;
end

function [t_values, X_values] = RK4_system(h, x0, y0, vx0, vy0, k1, k2, g)
    % Initialvärden
    t = 0;
    X = [x0; y0; vx0; vy0];
    
    % Lagra resultat
    t_values = t;
    X_values = X';
    
    while X(2) >= 0  % Fortsätt tills y < 0
        % Beräkna k-värden
        k1_vec = h * dynamics(X, k1, k2, g);
        k2_vec = h * dynamics(X + 0.5 * k1_vec, k1, k2, g);
        k3_vec = h * dynamics(X + 0.5 * k2_vec, k1, k2, g);
        k4_vec = h * dynamics(X + k3_vec, k1, k2, g);
        
        % Uppdatera tillståndet
        X = X + (k1_vec + 2 * k2_vec + 2 * k3_vec + k4_vec) / 6;
        t = t + h;
        
        % Lagra resultaten
        t_values = [t_values; t];
        X_values = [X_values; X'];
    end
end

function dXdt = dynamics(X, k1, k2, g)
    % Extrahera variabler
    vx = X(3);
    vy = X(4);
    speed = sqrt(vx^2 + vy^2);
    
    % Rörelseekvationer
    dxdt = vx;
    dydt = vy;
    dvxdt = -k1 * vx * speed;
    dvydt = -g - k2 * vy * speed;
    
    % Returnera derivator
    dXdt = [dxdt; dydt; dvxdt; dvydt];
end

function x_hit = interpolate_landing(X_values)
    % Hämta y-värden
    y_values = X_values(:,2);
    % Hitta index där y passerar genom noll
    idx = find(y_values <= 0, 1);
    
    % Säkerställ att vi har tillräckligt med punkter
    if idx < 4
        idx_range = 1:4;
    else
        idx_range = idx-3:idx;
    end
    
    % Extrahera punkter för interpolation
    x_pts = X_values(idx_range,1);
    y_pts = y_values(idx_range);
    
    % Utför kubisk interpolation
    p = polyfit(x_pts, y_pts, 3);
    % Hitta polynomets rötter
    roots_p = roots(p);
    % Välj den reella roten inom x-intervallet
    real_roots = roots_p(imag(roots_p) == 0);
    x_candidates = real_roots(real_roots >= min(x_pts) & real_roots <= max(x_pts));
    if isempty(x_candidates)
        x_hit = x_pts(end);
    else
        x_hit = x_candidates(1);
    end
end
