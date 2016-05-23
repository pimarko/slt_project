function [F] = frozen_bonds (curr_config, J, T, n, N)
    % Frozen bond matrix
    F = zeros(n, n);
    
    % Frozen bonds
    for i = 1:n
        for j = 1:n
            % Only both interact and if both are in the same state
            if is_neighbor(N, i, j) && J(i, j) > 0
                prob_frozen = 1 - exp((-J(i, j) / T) * (curr_config(i) == curr_config(j)));
                if rand() <= prob_frozen
                    F(i, j) = 1;
                end
            end
        end
    end
end

