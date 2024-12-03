function Main
    % Load model
    load('eiffel3.mat');  % Change this for other models
    n = length(xnod);  % Number of nodes
    
    % Initialize variables
    A_sparse = sparse(A);  % Convert matrix A to sparse
    
    % Perform sensitivity analysis and measure time
    methods = {'Naive', 'LU', 'Sparse Naive', 'Sparse LU'};
    results = zeros(4, 2);  % To store [min_index, max_index] for each method
    times = zeros(4, 1);    % To store execution time for each method
    
    % Method 1: Naive (dense matrix)
    tic;
    [min_node, max_node] = SensitivityAnalysis(A, 0);
    times(1) = toc;
    results(1, :) = [min_node, max_node];
    
    % Method 2: LU Decomposition (dense matrix)
    % [L, U]= lu(A) 
    tic;
    [L, U, P] = lu(A);  % Perform LU decomposition outside the loop
    [min_node, max_node] = SensitivityAnalysisLU(A, L, U, P);
    times(2) = toc;
    results(2, :) = [min_node, max_node];
    
    % Method 3: Sparse Naive
    tic;
    [min_node, max_node] = SensitivityAnalysis(A_sparse, 0);
    times(3) = toc;
    results(3, :) = [min_node, max_node];
    
    % Method 4: Sparse LU Decomposition
    tic;
    [L_sparse, U_sparse, P_sparse] = lu(A_sparse);  % LU for sparse matrix
    [min_node, max_node] = SensitivityAnalysisLU(A_sparse, L_sparse, U_sparse, P_sparse);
    times(4) = toc;
    results(4, :) = [min_node, max_node];
    
    % Display results
    disp('Results (Min, Max) for each method:');
    disp(table(methods', results(:,1), results(:,2), times, 'VariableNames', {'Method', 'MinNode', 'MaxNode', 'Time'}));

    % Save execution times to a text file
    fid = fopen('execution_times.txt', 'w');
    fprintf(fid, 'Method\t\tMinNode\tMaxNode\tTime (s)\n');
    for i = 1:4
        fprintf(fid, '%s\t%d\t%d\t%.6f\n', methods{i}, results(i,1), results(i,2), times(i));
    end
    fclose(fid);
    
    % =====================================
    % Plot the truss structure and mark the most/least sensitive nodes
    % =====================================
    figure;
    trussplot(xnod, ynod, bars);
    hold on;
    
    % Use results from method 1 (Naive method) for plotting
    most_sensitive_node = results(1, 2);  % Max sensitive node from Naive method
    least_sensitive_node = results(1, 1); % Min sensitive node from Naive method
    
    % Mark most sensitive node (red star) and least sensitive node (blue circle)
    plot(xnod(most_sensitive_node), ynod(most_sensitive_node), 'r*', 'MarkerSize', 10); % Mark most sensitive node
    plot(xnod(least_sensitive_node), ynod(least_sensitive_node), 'bo', 'MarkerSize', 10); % Mark least sensitive node
    
    title('Truss Structure with Sensitivity Markers');
    legend('Truss', 'Most Sensitive Node', 'Least Sensitive Node');
    hold off;
end

% Sensitivity Analysis for Naive or Sparse Naive methods
function [min_node, max_node] = SensitivityAnalysis(A, use_LU)
    n = size(A, 1) / 2;  % Number of nodes
    sensitivities = zeros(n, 1);  % Store sensitivities
    
    for j = 1:n
        % Create b vector for vertical force at node j
        b = zeros(2 * n, 1);
        b(2 * j) = -1;  % Apply vertical force
        
        % Solve the system
        xj = A \ b;  % Use direct solve for dense or sparse matrix
        
        % Compute the total displacement (Euclidean norm)
        sensitivities(j) = norm(xj);
    end
    
    % Find the most and least sensitive nodes
    [~, max_node] = max(sensitivities);
    [~, min_node] = min(sensitivities);
end

% Sensitivity Analysis using LU decomposition (dense or sparse)
function [min_node, max_node] = SensitivityAnalysisLU(A, L, U, P)
    n = size(A, 1) / 2;  % Number of nodes
    sensitivities = zeros(n, 1);  % Store sensitivities
    
    for j = 1:n
        % Create b vector for vertical force at node j
        b = zeros(2 * n, 1);
        b(2 * j) = -1;  % Apply vertical force
        
        % Solve the system using LU decomposition
        y = L \ (P * b); % L \ b
        xj = U \ y;
        
        % Compute the total displacement (Euclidean norm)
        sensitivities(j) = norm(xj);
    end
    
    % Find the most and least sensitive nodes
    [~, max_node] = max(sensitivities);
    [~, min_node] = min(sensitivities);
end

% Function to plot the truss structure
function A=trussplot(x,y,br,c)
% TRUSSPLOT Plots a truss
%
% TRUSSPLOT(X,Y,BARS) plots a truss with nodes in
% coordinates (X,Y) and bars between node indices
% given in BARS.
%
% The same plotstyles S as the PLOT command can be obtained
% with TRUSSPLOT(X,Y,BARS,S).
if (nargin < 4)
    c='k';
end
for k=1:length(br)
    plot(x(br(k,1:2)),y(br(k,1:2)),c); hold on
end
axis equal
hold off
end
