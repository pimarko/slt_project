function [J] = couplings (N, k, D, mean_D_normal, n)
    % Coupling matrix
    J = zeros(n, n);

    % Coupling for all neightbors of datapoints i
    parfor i = 1:n
        for j = 1:n
            if is_neighbor(N, i, j)
                J(i, j) = (1 / k) * exp(-D(i, j) / (2 * mean_D_normal^2));
            end
        end
    end
end

