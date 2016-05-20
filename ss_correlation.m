function [G] = swmc_sscorr (J, M, burns, q, n)
    % Configuration states: record inital spin configuration
	curr_config = sample_spins(q, n)';
    
    % Quantity to estimate
    C = zeros(n, n);
    
    % Swendsen-Wang Monte-Carlo estimator
    for iter = 1:M 
        % ---------------------------------------
        % 1) State transition: S_n -> S_n+1
        % ---------------------------------------
        [next_config, cluster_indices] = state_transition(curr_config, J, n);
        
        % ---------------------------------------
        % 2) Calc the quantitiy of interest
        % ---------------------------------------
        for i = 1:n
            for j = 1:n
                if cluster_indices(i) == cluster_indices(j) && iter > burns
                    C(i, j) = C(i, j) + 1;
                end
            end
        end
    end
    
    % Indicator function estimate
    C = C ./ (M - burns);

    % Spin-spin correlation function estimate
    G = (((q - 1) .* C) + 1) ./ q;
end

