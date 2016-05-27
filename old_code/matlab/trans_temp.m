function [T_init, T_final, T] = trans_temp (q, D, N, n, mean_D_normal)
    % Mean squared diff
    mean_D = mean(mean(D));
    
    T_trans = (exp(-mean_D / (2 * mean_D_normal^2))) / (4 * log10(1 + sqrt(q)));
    T_trans = 0.0346;
    %T_trans = exp(-1/2) / (4 * log10(1 + sqrt(q)));
    
    % Search space
    % Value provided by Viktor
    T_init = 4 * T_trans;
    T_final = 0.1 * T_trans;
    
    % Value provided by Marko
    %T_init = 1;
    %T_final = 0.001;
    
    % Iterating T
    T = T_init;
end