function [J] = couplings (X, N, k, D)
    % Number of datapoints
    n = size(X, 1);
    
    % Coupling matrix
    J = zeros(n, n);
    
    % D = Mean difference
    mean_diff = mean(mean(D));

    % Coupling for all neightbors of datapoints i
    for i = 1:n
        indices = N(N(:, i) > 0, i);
        for j = indices'
            J(i, j) = (1 / k) * exp(-D(i, j) / (2 * mean_diff^2));
        end
    end
end

