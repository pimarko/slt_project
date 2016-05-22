function [mean_k] = mean_neighbors(N)
    mean_k = 0;
    width = size(N, 2);
    for i = 1:width
        mean_k = mean_k + nnz(N(:, i));
    end
    mean_k = 2 * mean_k / width;
end

