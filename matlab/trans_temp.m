function [T_init, T_final, T] = trans_temp (q, D)
    % TODO: Extend to full equation ... ?
    T_trans = 1 / (4 * log(1 + sqrt(q))) * exp(-1/2);
    T_init = 4 * T_trans; % Value provided by Viktor
    T_final = 0.1 * T_trans; % Value provided by Viktor
    T = T_init;
end

