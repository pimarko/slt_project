function [J] = couplings (N, k, D, n)
    % Coupling matrix
    J = zeros(n, n);
    
    % D = Mean difference
    mean_diff = mean(mean(D));

    % Coupling for all neightbors of datapoints i
    for i = 1:n
        for j = 1:n
            if is_neighbor(N, i, j)
                J(i, j) = (1 / k) * exp(-D(i, j) / (2 * mean_diff^2));
            end
        end
    end
end

