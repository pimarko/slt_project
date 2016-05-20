function [next_config] = sample_cluster_comp (curr_config, nO_clusters, cluster_indices, q)
    next_config = curr_config;
    for i = 1:nO_clusters
        next_config(cluster_indices(cluster_indices == i)) = sample_spins(q, 1);
    end
end

