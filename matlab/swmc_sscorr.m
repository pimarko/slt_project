function [G] = swmc_sscorr (J, M, burns, q, n, T)
    % Configuration states: record inital spin configuration
    curr_config = sample_spins(q, n);

    % quantity to estimate
    C = zeros(n, n);
    
    % Swendsen-Wang Monte-Carlo estimator
    for iter = 1:M 
        % ---------------------------------------
        % 1) State transition: S_n -> S_n+1
        % ---------------------------------------
        [next_config, nO_clusters, cluster_indices] = state_transition(curr_config, J, T, n, q);
        
        % ---------------------------------------
        % 2) Calc the quantitiy of interest
        % ---------------------------------------
        % Calculate indicator function
        for i = 1:n
            for j = 1:n
                if cluster_indices(i) == cluster_indices(j) && iter > burns
                    C(i, j) = C(i, j) + 1;
                end
            end
        end
        
         % Switch the two spin configs
        curr_config = next_config;
    end
    
    % Indicator function estimate
    C = C ./ (M - burns);
    
    % spin-spin correlaction function
    G = (((q - 1) .* C) + 1) ./ q;
end

