function [F] = frozen_bonds (curr_config, J, T, n)
    % Frozen bond  matrix
    F = zeros(n, n);
    
    % Frozen bonds
    for i = 1:n
        for j = 1:n
            % Only between neighbors (J(i, j) > 0) and if both are in the same state
            if J(i, j) > 0 && curr_config(i) == curr_config(j)
                if rand() <= 1 - exp((-J(i, j) / T))
                    F(i, j) = 1;
                    F(j, i) = 1;
                end
            end
        end
    end
end

