function I = trapezoidal_2D(f, R, n)
    h = R / n;
    r_values = linspace(0, R, n + 1);
    g_values = f(r_values) .* r_values;
    I = (h / 2) * (g_values(1) + 2 * sum(g_values(2:end-1)) + g_values(end));
    I = 2 * pi * I; % Multiply by 2*pi for the polar integration
end
