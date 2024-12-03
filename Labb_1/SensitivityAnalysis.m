function [min_node, max_node] = SensitivityAnalysis(A, use_LU)
    % SensitivityAnalysis computes the sensitivity of each node to vertical forces
    % Inputs:
    % A - stiffness matrix (can be sparse)
    % use_LU - 0 for direct solve, 1 for LU decomposition
    % Outputs:
    % min_node - index of the least sensitive node
    % max_node - index of the most sensitive node
    
    n = size(A, 1) / 2;  % Number of nodes
    sensitivities = zeros(n, 1);  % Store sensitivities
    
    if use_LU
        [L, U, P] = lu(A);  % LU decomposition
    end
    
    for j = 1:n
        % Create b vector for vertical force at node j
        b = zeros(2 * n, 1);
        b(2 * j) = -1;  % Apply vertical force
        
        % Solve the system
        if use_LU
            y = L \ (P * b);
            xj = U \ y;
        else
            xj = A \ b;
        end
        
        % Compute the total displacement (Euclidean norm)
        sensitivities(j) = norm(xj);
    end
    
    % Find the most and least sensitive nodes
    [max_sensitivity, max_node] = max(sensitivities);
    [min_sensitivity, min_node] = min(sensitivities);
end
