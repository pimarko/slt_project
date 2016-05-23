function [neigbor] = is_neighbor (N, i, j)
    neigbor = false;
    if nnz(N(:, i) == j) > 0 && nnz(N(:, j) == i) > 0
        neigbor = true;
    end
end

