function x_impact = find_impact_point(t_vals, sol_vals, tolerance)
    % Extract x and y values
    x_vals = sol_vals(:, 1);
    y_vals = sol_vals(:, 2);
    
    % Find the last points where y is above and below zero
    idx = find(y_vals < 0, 1);  % First index where y < 0
    if isempty(idx) || idx < 2
        error('Impact point not found within simulation range or not enough points for interpolation');
    end
    
    % Select points around impact for interpolation (up to four)
    if idx >= 4
        x_interp = x_vals(idx-3:idx);
        y_interp = y_vals(idx-3:idx);
    else
        % Use fewer points if we have less than 4
        x_interp = x_vals(1:idx);
        y_interp = y_vals(1:idx);
    end
    
    % Perform cubic interpolation if we have 4 points, otherwise linear
    if length(x_interp) == 4
        p = polyfit(y_interp, x_interp, 3);  % Polynomial fit y -> x
        x_impact = polyval(p, 0);            % Evaluate at y = 0
    else
        % Linear interpolation as a fallback for fewer points
        x_impact = interp1(y_interp, x_interp, 0);
    end
    
    % Check error estimate
    if abs(x_vals(idx) - x_impact) > tolerance
        warning('Impact point may not meet tolerance requirements');
    end
end
