function [chi, noOfClusters] = swmc (J, X, s, M, burns, q, T)
    noOfDatapoints = size(X, 1);
    
    % Configuration matrix
    S = zeros(M + 1, noOfDatapoints);
    
    % record inital spin configuration
	S(1, :) = s';
        
    % Frozen matrix
    F = zeros(noOfDatapoints, noOfDatapoints);
    
    % estimated quantity m
    m = zeros(M, 1);
    %estimate_m_sqrd = zeros(M, 1);

    % SWMC
    for iter = 1:M 
        % Frozen bonds
        for i = 1:noOfDatapoints
            for j = 1:noOfDatapoints
                if J(i, j) > 0 && S(iter, i) == S(iter, j)
                    p_frozen = 1 - exp((-J(i, j) / T));
                    if rand() <= p_frozen
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
        
        % First need to calc n_max
        n_max = 0;
        for n = 1:q
            n_max = max(nnz(S(iter, :) == n), n_max);
        end
        
        % -------------
        % Calc quantities of interest
        % -------------
        m(iter, 1) = (q * n_max - noOfDatapoints) / ((q - 1) * noOfDatapoints);
        %estimate_m_sqrd(iter, 1) = m(iter, 1) * m(iter, 1);
    end
    
    % MC estimate <m> and <m>^2
    %estimate = sum(estimate_m(burns + 1:end)) / (M - burns);
    %estimate_m_sqrd = sum(estimate_m_sqrd(burns + 1:end)) / (M - burns);

    % Calculate Chi
    %chi = (noOfDatapoints / T) * (estimate_m_sqrd - estimate * estimate);
    chi = (noOfDatapoints / T) * var(m(burns + 1:end));
 end

