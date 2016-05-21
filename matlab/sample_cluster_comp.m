function [next_config] = sample_cluster_comp (curr_config, nO_clusters, cluster_indices, q)
    next_config = zeros(1, size(curr_config, 2));
    for i = 1:nO_clusters
        next_config(cluster_indices == i) = sample_spins(q, 1);
    end
end

