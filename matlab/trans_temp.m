function [T_init, T_final, T] = trans_temp (q, D, N, n)
    % Mean distance between al neighboring pairs
%     count = 1;
%     mean_D = 0;
%     for i = 1:n
%         indices = N(N(:, i) > 0, i);
%         for j = indices'
%             mean_D = mean_D + D(i, j);
%             count = count + 1;
%         end
%     end
%     
%     mean_D = mean_D / count;

    % Estimated transition temparature
%    T_trans = 1 / (4 * log(1 + sqrt(q))) * exp(-(mean_D) / (2 * mean(mean(D))));
    T_trans = 1 / (4 * log(1 + sqrt(q))) * exp(-1/2);
    
    % Search space
    T_init = 4 * T_trans; % Value provided by Viktor
    T_final = 0.1 * T_trans; % Value provided by Viktor
    T = T_init;
end

