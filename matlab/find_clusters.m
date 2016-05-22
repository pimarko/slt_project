function [nO_Clusters, cluster_indices] = find_clusters (G, n)
    % Clusters linkage matrix
    cluster_linkage = eye(n, n);
    
    % Find linked vertices
    for i = 1:n
        for j = 1:n
            if G(i, j) > 0.5
                cluster_linkage(i, j) = 1;
            end
        end
    end
    
    % Capture peripherie
    for i = 1:n
        p_max_neighbor = find(G(i, :) == max(G(i, :)), 1, 'first');
        cluster_linkage(i, p_max_neighbor) = 1;
    end
    
    % Find cluster_linkage
	[nO_Clusters, cluster_indices] = graphconncomp(sparse(cluster_linkage), 'weak', true);
end

