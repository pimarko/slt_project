addpath('nifti');
clearvars -except data bvecs

q       = 10; % Spin states
M       = 100; % Monte carlo samples to draw
burns   = 10; % Monte carlo burn samples

i_init  = 135;
j_init  = 76;
k_init  = 74;

i_max   = 21; % Patch size in i
j_max   = 21; % Patch size in j
k_max   = 21; % Patch size in k

n       = i_max * j_max * k_max; % Number of datapoints
k       = 26; % Number of nearest neighbors
eta     = 0.95; % Exponential cooling

% Read data (only of not already read)
if any(strcmp(who, 'data')) == 0
    data = load_nii('data/diff_data.nii.gz');
    bvecs = dlmread('data/bvecs.txt');
end

% -----------------------------------------------------------------
% 1) Preprocessing: Load data and build neighborhood matrix
% -----------------------------------------------------------------
[X, N, D, mean_D_normal, coordinate_map] = read_data(data, bvecs, k, i_init, i_max, j_init, j_max, k_init, k_max);

% -----------------------------------------------------------------
% 2) Calculate couplings between neighbors
% -----------------------------------------------------------------
J = couplings(N, mean_neighbors(N), D, mean_D_normal, n);

% Temperature estimate for ferro -> para
[T_init, T_final, T] = trans_temp(q, D, N, n, mean_D_normal);

% Variables to keep track of the values
Ts = []; chis = []; clusters = []; ms = [];

% -----------------------------------------------------------------
% 3) Locate the super-paramagnetic phase
% -----------------------------------------------------------------
% iter = 0;
% while T > T_final
%     % SWMC for Chi
%     [chi, m, nO_clusters] = swmc_chi(J, M, N, burns, q, T, n);
%     
%     % Track values
%     chis = [chi; chis]; clusters = [nO_clusters; clusters]; Ts = [T; Ts]; ms = [m; ms];
%     
%     % Exponential cooling
%     T = T_init * (eta ^ iter);
%     iter = iter + 1
% end

% Locate temperature in the paramagnetic region
T = Ts(find(chis == max(chis), 1, 'first') + 1)
T = 0.15;

% -----------------------------------------------------------------
% 3) Calc the spin-spin correlation matrix
% -----------------------------------------------------------------
G = swmc_sscorr(N, J, M, burns, q, n, T);

% -----------------------------------------------------------------
% 4) Find clusters
% -----------------------------------------------------------------
[nO_clusters, cluster_indices] = find_clusters(G, N, n);

% -----------------------------------------------------------------
% 5) Plot
% -----------------------------------------------------------------
scatter3(coordinate_map(:, 1), coordinate_map(:, 2), coordinate_map(:, 3), 800, cluster_indices, 'filled', 'square');

% Chi
fig=figure; 
hax=axes;  
hold on 
plot(Ts, chis) 
line([T T],get(hax,'YLim'),'Color',[1 0 0])
xlabel('T');
ylabel('Chi');

% <m>
figure;
plot(Ts, ms)
xlabel('T');
ylabel('<m>');

% G
figure;
hist(G(G >= 1/q))
xlabel('G');

%figure;
%plot(Ts, clusters);