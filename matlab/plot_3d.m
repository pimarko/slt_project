voxel_count = zeros(nO_clusters, 1);

% Count voxels per clusters
for i = 1:n
   voxel_count(cluster_indices(i)) = voxel_count(cluster_indices(i)) + 1;
end

% New coordinate map
clrd_coordinates = [];
clrd_clusters = [];

for i = 1:n
    cluster_index = cluster_indices(i);
    % Only clusters which ahve more than x voxels
    if voxel_count(cluster_index) > 180
        clrd_coordinates = [clrd_coordinates; coordinate_map(i, :)];
        clrd_clusters = [clrd_clusters, cluster_index];
    end
end

% 3d scatter plot
scatter3(clrd_coordinates(:, 1), clrd_coordinates(:, 2), clrd_coordinates(:, 3), 400, clrd_clusters, 'filled', 'square');