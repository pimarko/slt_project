function [J] = couplings (X, N, k, D)
    % Coupling matrix
    noOfDatapoints = size(X, 1);
    J = zeros(noOfDatapoints, noOfDatapoints);
    
    % D = Mean difference
    mean_diff = mean(mean(D));
    
    %N_n = size(N, 2);
    for i = 1:noOfDatapoints
        indices = N(N(:, i) > 0, i);
        for j = indices'
            J(i, j) = (1 / k) * exp(-D(i, j) / (2 * mean_diff^2));
        end
    end
end

