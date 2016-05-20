function [next_config, nO_clusters, cluster_indices] = state_transition (curr_config, J, T, n, q)
    % Calculate Frozen bonds
    F = frozen_bonds(curr_config, J, T, n);

    % Build graph and locate connected components
    [nO_clusters, cluster_indices] = graphconncomp(sparse(F), 'weak', true);

    % Assign a new spin value to each vertice in the connected component
    next_config = sample_cluster_comp(curr_config, nO_clusters, cluster_indices, q);
end

