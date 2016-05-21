function [chi, nO_clusters] = swmc_chi (J, M, burns, q, T, n)
    % Configuration states: record inital spin configuration
    curr_config = sample_spins(q, n);

    % Quantity to estimate
    m = zeros(M, 1);

    % Swendsen-Wang Monte-Carlo estimator
    for iter = 1:M
        % ---------------------------------------
        % 1) State transition: S_n -> S_n+1
        % ---------------------------------------
        [next_config, nO_clusters, cluster_indices] = state_transition(curr_config, J, T, n, q);
        
        % ---------------------------------------
        % 2) Calc the quantitiy of interest
        % ---------------------------------------
        % First need to calc n_max
        n_max = 0;
        for i = 1:q
            n_max = max(nnz(curr_config(:) == i), n_max);
        end
        
        % Now calc m
        m(iter) = (q * n_max - n) / ((q - 1) * n);
        
        % Switch the two spin configs
        curr_config = next_config;
    end

    % Estimate X via SWMC <X>
    if burns == 0
        burns = 1;
    end
    
    chi = (n / T) * var(m(burns:end));  
 end

