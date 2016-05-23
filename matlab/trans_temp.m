function [T_init, T_final, T] = trans_temp (q, D, N, n)
    T_trans = exp(-1/2) / (4 * log(1 + sqrt(q)));
    
    % Search space
    % Value provided by Viktor
    %T_init = 4 * T_trans;
    %T_final = 0.1 * T_trans;
    
    % Value provided by Marko
    T_init = T_trans + 1;
    T_final = T_trans - T_trans / 5;
    
    % Iterating T
    T = T_init;
end