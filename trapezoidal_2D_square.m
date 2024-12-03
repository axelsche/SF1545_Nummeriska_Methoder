function I = trapezoidal_2D_square(f, L, n)
    h = L / n;
    x = linspace(-L/2, L/2, n + 1);
    y = linspace(-L/2, L/2, n + 1);
    
    % Initialize sum with edge cases
    I = 0;
    for i = 1:n+1
        for j = 1:n+1
            weight = 1;
            if i == 1 || i == n+1, weight = weight / 2; end
            if j == 1 || j == n+1, weight = weight / 2; end
            I = I + weight * f(x(i), y(j));
        end
    end
    
    % Multiply by area element (h^2) to get the final result
    I = h^2 * I;
end
