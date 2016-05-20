addpath('nifti');

q = 10; % Spin states
M = 100; % Monte carlo samples to draw
burns =1; % Monte carlo burn samples
i_max = 9; % Patch size in i
j_max = 9; % Patch size in j
k_max = 9; % Patch size in k
n = i_max * j_max * k_max; % Number of datapoints
k = 26; % Number of nearest neighbors
eta = 0.95; % Exponential cooling

% Read data (only of not already read)
if any(strcmp(who, 'data')) == 0
    data = load_nii('data/diff_data.nii.gz');
    bvecs = dlmread('data/bvecs.txt');
end

% -----------------------------------------------------------------
% 1) Preprocessing: Load data and build neighborhood matrix
% -----------------------------------------------------------------
[X, N, D, coordinate_map] = read_data(data, bvecs, k, i_max, j_max, k_max);

% -----------------------------------------------------------------
% 2) Calculate couplings between neighbors
% -----------------------------------------------------------------
J = couplings(X, N, mean_neighbors(N), D);

% Temperature estimate for ferro -> para
[T_init, T_final, T] = trans_temp(q, D);

% Variables to keep track of the values
Ts = []; chis = []; clusters = [];

% -----------------------------------------------------------------
% 3) Locate the super-paramagnetic phase
% -----------------------------------------------------------------
iter = 0;
while T > T_final
    % SWMC for Chi
    [chi, nO_clusters] = swmc_chi(J, X, M, burns, q, T, n);
    
    % Track values
    chis = [chi; chis]; clusters = [nO_clusters; clusters]; Ts = [T;Ts];
    
    % Exponential cooling
    T = T_init * (eta ^ iter);
    iter = iter + 1;
end

% Locate temperature in the paramagnetic region
T = Ts(find(chis == max(chis), 1, 'first'));
T = 0.08
% -----------------------------------------------------------------
% 3) Calc the spin-spin correlation matrix
% -----------------------------------------------------------------
G = swmc_sscorr(J, M, burns, q, n, T);

% -----------------------------------------------------------------
% 4) Find clusters
% -----------------------------------------------------------------
[noOfClusters, cluster_indices] = find_clusters(G, n);

% -----------------------------------------------------------------
% 5) Plot
% -----------------------------------------------------------------
scatter3(coordinate_map(:, 1), coordinate_map(:, 2), coordinate_map(:, 3), 800, cluster_indices, 'filled', 'square');

fig=figure; 
hax=axes;  
hold on 
plot(Ts, chis) 
line([T T],get(hax,'YLim'),'Color',[1 0 0])

figure;
plot(Ts, clusters);