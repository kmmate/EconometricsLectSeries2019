function [normdiff_normal_leverage, normdiff_normal_uniform, ...
          normdiff_skewed_leverage, normdiff_skewed_uniform] = ...
          q2c_p2_simulation(X_normal,...
                            X_skewed,...
                            m,...
                            reps)
%% 

%% Moments of Xs
SX_normal = X_normal' * X_normal;
SX_skewed = X_skewed' * X_skewed;
%% Distributions used for resampling
% computing leverage scores
[n, k] = size(X_normal);
h_normal = efficient_H_diagonal(X_normal);
h_skewed = efficient_H_diagonal(X_skewed);
p_normal = h_normal / k;
p_skewed = h_skewed / k;
pd_leverage_normal = makedist('Multinomial', 'probabilities', p_normal);
pd_leverage_skewed = makedist('Multinomial', 'probabilities', p_skewed);
pd_uniform = makedist('Multinomial', 'probabilities', ones(n,1)/n);

%% Simulation
normdiff_normal_leverage = zeros(reps, 1);
normdiff_normal_uniform = zeros(reps, 1);
normdiff_skewed_leverage = zeros(reps, 1);
normdiff_skewed_uniform = zeros(reps, 1);
for rep = 1:reps
    % --- Normal
    
    % - resample with leverage
    leverage_idx = random(pd_leverage_normal, m, 1);
    Xnormal_leverage = X_normal(leverage_idx, :);
    % norm of difference
    SXnormal_leverage = Xnormal_leverage' * Xnormal_leverage;
    normdiff_normal_leverage(rep) = norm(SXnormal_leverage - SX_normal, 'fro');
    
    % - resample with uniform
    uniform_idx = random(pd_uniform, m, 1);
    Xnormal_uniform = X_normal(uniform_idx, :);
    % norm of difference
    SXnormal_uniform = Xnormal_uniform' * Xnormal_uniform;
    normdiff_normal_uniform(rep) = norm(SXnormal_uniform - SX_normal, 'fro');
    
    
    % --- Skewed
    
    % - resample with leverage
    leverage_idx = random(pd_leverage_skewed, m, 1);
    Xskewed_leverage = X_skewed(leverage_idx, :);
    % norm of difference
    SXskewed_leverage = Xskewed_leverage' * Xskewed_leverage;
    normdiff_skewed_leverage(rep) = norm(SXskewed_leverage - SX_skewed, 'fro');
    
    % - resample with uniform
    uniform_idx = random(pd_uniform, m, 1);
    Xskewed_uniform = X_skewed(uniform_idx, :);
    % norm of difference
    SXskewed_uniform = Xskewed_uniform' * Xskewed_uniform;
    normdiff_skewed_uniform(rep) = norm(SXskewed_uniform - SX_skewed, 'fro');    
end