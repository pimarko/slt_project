function [X, N, N_indices, D] = read_data (data, noOfNeighbors, i_max, j_max, k_max)
    bvecs = dlmread('data/bvecs.txt');
    
    % Count non-zero rows
    noOfNZRows = size(bvecs(any(bvecs, 2), :), 1);
    
    % Data Matrix: Preallocation
    noOfDatapoints = i_max * j_max * k_max;
    X = zeros(noOfDatapoints, noOfNZRows);
    N_indices = zeros(noOfDatapoints, 3);
    
    % Neigbor mapping
    N_map = zeros(i_max, j_max, k_max);
    
    % counts the number of samples
    index = 1;
    
    % Extract voxel with direction vector n
    for i = 1:i_max 
        for j = 1:j_max
            for k = 1:k_max
                index_n = 1;
                for n = 1:164
                    % Only non-zero directions
                    if bvecs(n, 1) ~= 0 && bvecs(n, 2) ~= 0 && bvecs(n, 3) ~= 0
                        % Signal strength in direction n
                        X(index, index_n) = data.img(i, j, k, n);
                        index_n = index_n + 1;
                    end
                end
                
                % Reference neighbors
                N_map(i, j, k) = index;
                N_indices(index, : ) = [i, j ,k];
                
                index = index + 1;
            end
        end
    end

    
    % Build kd-tree to search for k neighbors
    Mdl = KDTreeSearcher(N_indices);
    
    % Neighbor and distance matrix
    N = zeros(noOfNeighbors, noOfDatapoints);
    D = squareform(pdist(N_indices), 'tomatrix');
    
    for i = 1:i_max 
        for j = 1:j_max
            for k = 1:k_max
                indices = knnsearch(Mdl, [i, j, k], 'K', noOfNeighbors);
                N(:, N_map(i, j, k)) = indices;
            end
        end
    end
    
    % Restrict to mutual neighbors only
    for n = 1:noOfDatapoints
        for j = 1:noOfNeighbors
            if N(j, n) > 0 && nnz(N(:, N(j, n)) == n) == 0
                N(j, n) = 0;
            end
        end
    end
end

