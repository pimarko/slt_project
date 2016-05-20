addpath('nifti');

q = 10; % Spin states
M = 100; % Monte carlo samples to draw
burns = 10; % Monte carlo burn samples
i_max = 5; % Patch size in i
j_max = 5; % Patch size in j
k_max = 5; % Patch size in k
k = 5; % Number of nearest neighbors
eta = 0.97;

% Read data (only of not already read)
if any(strcmp(who, 'data')) == 0
    data = load_nii('data/diff_data.nii.gz');
end

% -------------
% Preprocessing
% -------------
[X, N, N_indices, D] = read_data(data, k, i_max, j_max, k_max);

% -------------
% Couplings
% -------------
% k = mean no. of neighbors
mean_k = 0;
width = size(N, 2);
for i = 1:width
    mean_k = mean_k + nnz(N(:, i));
end
mean_k = mean_k / width;
J = couplings(X, N, mean_k, D);

% -------------
% Inital spin configuration
% -------------
s = randi([1, q], size(X, 1), 1);

% -------------
% Temperature estimates
% -------------
T_estimate = 1 / (4 * log(1 + sqrt(q))) * exp(-1/2);
T_inital = T_estimate + 1; % Werte von marko
T_final = T_estimate - (T_estimate / 2); % Werte von marko
T = T_inital;

% Vars to keep track of the values
Ts = [];
chis = [];
cluster_trace = [];

% Simluate temperature
iter = 0;
while T > T_final
    % -------------
    % SWMC for Chi
    % -------------
    [chi, noOfClusters] = swmc(J, X, s, M, burns, q, T);
    chis = [chis;chi];
    cluster_trace = [cluster_trace; noOfClusters];
    
    % Update temperature
    Ts = [Ts;T];
    
    % Exponential cooling
    T = T_inital * (eta ^ iter);
    iter = iter + 1;
end

% Temperature in the superpara. region
T = Ts(find(chis == max(chis), 1, 'first'));
%T = 0.2;

% -------------
% Spin-spin correlation
% -------------
G = ss_correlation(J, T, M, burns, q, size(X, 1));

% -------------
% Find clusters
% -------------
[noOfClusters, cluster_indices] = find_clusters(G, N, q);

% -------------
% Plot
% -------------
scatter3(N_indices(:, 1), N_indices(:, 2), N_indices(:, 3), 30, cluster_indices);

%fig=figure; 
%hax=axes;  
%hold on 
%plot(Ts, chis) 
%line([T T],get(hax,'YLim'),'Color',[1 0 0])

%figure;
%plot(Ts, cluster_trace);