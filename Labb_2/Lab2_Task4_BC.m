% varpa_delb.m
% Beräkning av kastvinklar som träffar pinnen vid x = 20 meter
% och plottar både höga och låga banor.

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
h0 = 0.5;               % Start höjd [m]

% Toleranser
tolerance_alpha = 1e-6; % Tolerans för kastvinkeln
tolerance_x = 1e-6;     % Tolerans för nedslagsplatsen
h_initial = 0.01;       % Initial steglängd för RK4

% Målvärde för x
x_target = 20;          % Pinnen finns vid x = 20 meter

% Initiala gissningar för låga banan
alpha_low1 = 0.1;       % Lägre gissning för låg bana [rad]
alpha_low2 = 0.6;       % Högre gissning för låg bana [rad]

% Initiala gissningar för höga banan
alpha_high1 = 0.8;      % Lägre gissning för hög bana [rad]
alpha_high2 = 1.2;      % Högre gissning för hög bana [rad]

% Beräkna låga kastvinkeln med intervallhalvering
alpha_low = bisection_method(@(alpha) trajectory_function(alpha, v0, h0, x_target, m, g, k1, k2, h_initial, tolerance_x), alpha_low1, alpha_low2, tolerance_alpha);

% Beräkna höga kastvinkeln med intervallhalvering
alpha_high = bisection_method(@(alpha) trajectory_function(alpha, v0, h0, x_target, m, g, k1, k2, h_initial, tolerance_x), alpha_high1, alpha_high2, tolerance_alpha);

% Beräkna banorna för de funna vinklarna
[~, ~, X_low] = calculate_trajectory(m, g, k1, k2, v0, alpha_low, h0, h_initial, tolerance_x);
[~, ~, X_high] = calculate_trajectory(m, g, k1, k2, v0, alpha_high, h0, h_initial, tolerance_x);

% Visa resultaten
fprintf('Låg kastvinkel alpha_low: %.8f rad (%.6f grader)\n', alpha_low, rad2deg(alpha_low));
fprintf('Hög kastvinkel alpha_high: %.8f rad (%.6f grader)\n', alpha_high, rad2deg(alpha_high));

% Plotta banorna
figure;
plot(X_low(:,1), X_low(:,2), 'b-', 'LineWidth', 1.5);
hold on;
plot(X_high(:,1), X_high(:,2), 'r-', 'LineWidth', 1.5);
plot(x_target, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');  % Pinnen vid x = 20
xlabel('x [m]');
ylabel('y [m]');
title('Varpastenens banor som träffar pinnen vid x = 20 m');
legend('Låg bana', 'Hög bana', 'Pinnen', 'Location', 'Best');
grid on;

%% Osäkerhetsberäkningar

% Grundvärden
alpha_low_base = alpha_low;
alpha_high_base = alpha_high;
v0_base = v0;

% Omvandla grundvinklar till grader
alpha_low_deg = rad2deg(alpha_low_base);
alpha_high_deg = rad2deg(alpha_high_base);

% Definiera ±5% variationer i grader
delta_alpha_low_deg = 0.05 * alpha_low_deg;
delta_alpha_high_deg = 0.05 * alpha_high_deg;

% Skapa variationsvektorer i grader
alpha_variations_low_deg = [alpha_low_deg - delta_alpha_low_deg, alpha_low_deg + delta_alpha_low_deg];
alpha_variations_high_deg = [alpha_high_deg - delta_alpha_high_deg, alpha_high_deg + delta_alpha_high_deg];

% Omvandla tillbaka till radianer
alpha_variations_low = deg2rad(alpha_variations_low_deg);
alpha_variations_high = deg2rad(alpha_variations_high_deg);

% Definiera ±5% variationer för v0
delta_v0 = 0.05 * v0_base;
v0_variations = [v0_base - delta_v0, v0_base + delta_v0];

% Initialisera arrayer för att lagra nedslagsplatser
x_hits_low = [];
x_hits_high = [];

% Beräkna x_hit för låg bana (alla kombinationer av störda variabler)
for alpha_var = alpha_variations_low
    for v0_var = v0_variations
        f = trajectory_function(alpha_var, v0_var, h0, x_target, m, g, k1, k2, h_initial, tolerance_x);
        x_hit = f + x_target;
        x_hits_low = [x_hits_low; x_hit];
    end
end

% Beräkna osäkerhet för låg bana
uncertainty_low = (max(x_hits_low) - min(x_hits_low)) / 2;

% Beräkna x_hit för hög bana (alla kombinationer av störda variabler)
for alpha_var = alpha_variations_high
    for v0_var = v0_variations
        f = trajectory_function(alpha_var, v0_var, h0, x_target, m, g, k1, k2, h_initial, tolerance_x);
        x_hit = f + x_target;
        x_hits_high = [x_hits_high; x_hit];
    end
end

% Beräkna osäkerhet för hög bana
uncertainty_high = (max(x_hits_high) - min(x_hits_high)) / 2;

% Visa osäkerheterna
fprintf('Osäkerhet i nedslagsplats för låg kastvinkel: %.6f meter\n', uncertainty_low);
fprintf('Osäkerhet i nedslagsplats för hög kastvinkel: %.6f meter\n', uncertainty_high);

% Avgör vilket kast som är säkrare
if uncertainty_low < uncertainty_high
    fprintf('Det låga kastet är säkrare.\n');
elseif uncertainty_low > uncertainty_high
    fprintf('Det höga kastet är säkrare.\n');
else
    fprintf('Båda kasten har samma osäkerhet.\n');
end

% Funktioner
function f = trajectory_function(alpha, v0, h0, x_target, m, g, k1, k2, h_initial, tolerance_x)
    % Beräknar differensen mellan nedslagsplatsen och målet x_target
    [x_hit, ~, ~] = calculate_trajectory(m, g, k1, k2, v0, alpha, h0, h_initial, tolerance_x);
    f = x_hit - x_target;
end

function alpha_root = bisection_method(f, alpha_left, alpha_right, tol)
    % Intervallhalveringsmetoden för att hitta roten till f(alpha) = 0
    f_left = f(alpha_left);
    f_right = f(alpha_right);
    
    if f_left * f_right > 0
        error('Funktionen har samma tecken vid intervallets ändpunkter.');
    end
    
    while abs(alpha_right - alpha_left) > tol
        alpha_mid = (alpha_left + alpha_right) / 2;
        f_mid = f(alpha_mid);
        
        if f_left * f_mid < 0
            alpha_right = alpha_mid;
            f_right = f_mid;
        else
            alpha_left = alpha_mid;
            f_left = f_mid;
        end
    end
    
    alpha_root = (alpha_left + alpha_right) / 2;
end

function [x_hit, error_estimate, X1] = calculate_trajectory(m, g, k1, k2, v0, alpha, h0, h_initial, tolerance_x)
    % Beräkna banan för en given vinkel alpha
    % Begynnelsevärden
    vx0 = v0 * cos(alpha);
    vy0 = v0 * sin(alpha);
    x0 = 0;
    y0 = h0;
    
    % Initial steglängd
    h_current = h_initial;
    % Initialisera feluppskattning
    error_estimate = Inf;
    
    while error_estimate > tolerance_x
        % Lös med steglängd h_current
        [t1, X1] = RK4_system(h_current, x0, y0, vx0, vy0, k1, k2, g);
        x1_hit = interpolate_landing(X1);
        
        % Lös med dubbel steglängd
        [t2, X2] = RK4_system(2*h_current, x0, y0, vx0, vy0, k1, k2, g);
        x2_hit = interpolate_landing(X2);
        
        % Beräkna feluppskattning
        error_estimate = abs(x1_hit - x2_hit);
        
        if error_estimate > tolerance_x
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
