function [F] = frozen_bonds (curr_config, J, T, n, N)
    % Frozen bond  matrix
    F = zeros(n, n);
    
    % Frozen bonds
    for i = 1:n
        indices = N(N(:, i) > 0, i);
        for j = indices'
            % Only if both are in the same state
            if rand() <= (1 - exp((-J(i, j) / T) * (curr_config(i) == curr_config(j))))
                F(i, j) = 1;
                F(j, i) = 1;
            end
        end
    end
end

