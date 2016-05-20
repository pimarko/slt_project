function [noOfClusters, cluster_indices] = find_clusters (G, N, q)
    % Clusters
    noOfDatapoints = size(N, 2);
    clusters = zeros(noOfDatapoints, noOfDatapoints);
    
    % Find linked vertices
    for i = 1:noOfDatapoints
        for j = 1:noOfDatapoints
            if G(i, j) > 0.5
                clusters(i, j) = 1;
                clusters(j, i) = 1;
            end
        end
    end
    
    % Capture peripherie
    for i = 1:noOfDatapoints
        p_max_neighbor = find(G(i, :) == max(G(i, :)), 1, 'first');
        clusters(i, p_max_neighbor) = 1;
        clusters(p_max_neighbor, 1) = 1;
    end
    
    % Find clusters
	[noOfClusters, cluster_indices] = graphconncomp(sparse(clusters));
end

