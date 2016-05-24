function [X, N, D, index_map] = read_data (data, bvecs, nO_neighbors, i_max, j_max, k_max)
    % Count non-zero rows
    nO_nz_rows = size(bvecs(any(bvecs, 2), :), 1);
    
    % Number of data points
    n = i_max * j_max * k_max;
    
    % Data Matrix: Preallocation
    X = zeros(n, nO_nz_rows);
    
    % Index to coordinates (i, j, k) mapping
    index_map = zeros(n, 3);
    
    % Coordinates (i, j, k) to index mapping
    coordinate_map = zeros(i_max, j_max, k_max);
    
    % Running datapoint v_index
    index = 1;
    
    % Extract voxel with direction vector n
    for i = 1:i_max 
        for j = 1:j_max
            for k = 1:k_max
                index_direction = 1;
                for d = 1:164
                    % Only non-zero directions
                    if bvecs(d, 1) ~= 0 || bvecs(d, 2) ~= 0 || bvecs(d, 3) ~= 0
                        % Signal strength in direction n
                        X(index, index_direction) = data.img(i, j, k, d);
                        index_direction = index_direction + 1;
                        %sprintf('%i, %i, %i, %i', i, j, k, d)
                    end
                end
                
                % Update coordinate / index map
                coordinate_map(i, j, k) = index;
                index_map(index, : ) = [i, j ,k];
                
                index = index + 1;
            end
        end
    end

    % Build kd-tree to search for k neighbors
    Mdl = KDTreeSearcher(index_map);
    
    % Neighbor and distance matrix
    N = zeros(nO_neighbors, n);
    D = squareform(pdist(X), 'tomatrix');
    D = D.^2;
    
    % Query k-nearest neighbors
    for i = 1:i_max 
        for j = 1:j_max
            for k = 1:k_max
                indices = knnsearch(Mdl, [i, j, k], 'K', nO_neighbors);
                N(:, coordinate_map(i, j, k)) = indices;
            end
        end
    end
    
    % Restrict to mutual neighbors only
    for i = 1:n
        for j = 1:nO_neighbors
            if N(j, i) > 0 && nnz(N(:, N(j, i)) == i) == 0
                N(j, i) = 0;
            end
        end
    end
end

