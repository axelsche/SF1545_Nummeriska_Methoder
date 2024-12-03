function I = simpson_1D(f, a, b, n)
    if mod(n, 2) ~= 0
        error('n must be even for Simpson''s rule');
    end
    h = (b - a) / n;
    x = linspace(a, b, n + 1);
    y = f(x);
    I = (h / 3) * (y(1) + 4 * sum(y(2:2:end-1)) + 2 * sum(y(3:2:end-2)) + y(end));
end
