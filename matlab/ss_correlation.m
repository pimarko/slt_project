function [G] = ss_correlation (J, T, M, burns, q, noOfDatapoints)

    % Configuration matrix
    S = zeros(M + 1, noOfDatapoints);
    
    % record inital spin configuration
	S(1, :) = randi([1, q], 1, 1);
    
    % Frozen matrix
    F = zeros(noOfDatapoints, noOfDatapoints);
    
    % quantity to estimate
    C = zeros(noOfDatapoints, noOfDatapoints);
    
    % SWMC
    for iter = 1:M 
        % Frozen bonds
        for i = 1:noOfDatapoints
            for j = 1:noOfDatapoints
                if J(i, j) > 0 && S(iter, i) == S(iter, j)
                    p_frozen = exp((-J(i, j) / T));
                    if rand() > p_frozen
                        F(i, j) = 1;
                        F(j, i) = 1;
                    end
                end
            end
        end
        
        % Build graph and locate connected components
        [noOfClusters, cluster_indices] = graphconncomp(sparse(F));
        
        % Asign new values to all vertices in connected components
        for i = 1:noOfClusters
            s = randi([1, q], 1, 1); % 1 new spin state
            S(iter + 1, cluster_indices(cluster_indices == i)) = s; % trace
        end
        
        % Calculate indicator function
        for i = 1:noOfDatapoints
            for j = 1:noOfDatapoints
                if cluster_indices(i) == cluster_indices(j) && iter > burns
                    C(i, j) = C(i, j) + 1;
                end
            end
        end
    end
    
    % Indicator function estimate
    C = C ./ (M - burns);
    
    % spin-spin correlaction function
    G = (((q - 1) .* C) + 1) ./ q;
end

