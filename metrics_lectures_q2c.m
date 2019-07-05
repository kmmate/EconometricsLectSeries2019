%% Tinbergen Econometrics Lecture Series 2019, Exam, Question 2 / c
%% C (i) Computing p
fprintf('\n ==============     C (i)     =========================\n')
% setup
rng(20190704)  % set random seed
n = 10^6;
k = 10;
% draw data
X = randn(n, k);  % data from iid N(0,1)

% hat matrix diagonal elements
h = efficient_H_diagonal(X);
p = h / k;

% test uniformity
pd = makedist('Multinomial', 'probabilities', p);
idx_sample = random(pd, 1 * n, 1);
test_pd = makedist('Uniform', 'Lower', 1, 'Upper', n);
[hypothesis, p_value, test_stat] = kstest(idx_sample, 'CDF', test_pd);
fprintf('- Uniformity test p-value: %.4f', p_value)
%% C (ii) Approximate p with p_tilde
fprintf('\n ==============     C (ii)     =========================')
r = 4;
core = (n ^ (2/r)) * (k^(1+2/r)) * log(k);
m_range = round(core * [1 1.5 2]);%round(core * [0.01 0.1 1 2]);
for m=m_range
   % computing p_tilde
   Xtilde = uniform_sampling(X, m); % subsample X
   h_tilde = efficient_H_diagonal(Xtilde);
   p_tilde = h_tilde / k;
   
   % test uniformity
   pd = makedist('Multinomial', 'probabilities', p_tilde);
   idx_sample = random(pd, m, 1);
   figure
   cdfplot(idx_sample)
   test_pd = makedist('Uniform', 'Lower', 1, 'Upper', m);
   [hypothesis, p_value, test_stat] = kstest(idx_sample, 'CDF', test_pd); 
   fprintf('\n --- m = %d', m)   
   fprintf('- \n Uniformity test p-value: %.4f', p_value)
end

%% C (iii) Skewed X
% fprintf('\n ==============     C (iii)     =========================')
% pdX = makedist('Gamma', 'a', 5, 'b', 14);
% sum_pvalue_skewed = 0;
% sum_pvalue_tilde = 0;
% for rep = 1:50
% fprintf('\n Rep = %d', rep)
% X_skewed = random(pdX, n, k);
% 
% % --- True p
% % hat matrix diagonal elements
% h = efficient_H_diagonal(X_skewed);
% p_skewed = h / k;
% % test uniformity
% pd = makedist('Multinomial', 'probabilities', p_skewed);
% idx_sample = random(pd, 1 * n, 1);
% test_pd = makedist('Uniform', 'Lower', 1, 'Upper', n);
% [hypothesis, p_value, test_stat] = kstest(idx_sample, 'CDF', test_pd);
% fprintf('\n ---- Original p_skewed')
% fprintf('\n- Uniformity test p-value: %.4f', p_value)
% sum_pvalue_skewed = sum_pvalue_skewed + p_value;
% % --- Approximated p
% fprintf('\n---- Approximated p_skewed')
% r = 4;
% core = (n ^ (2/r)) * (k^(1+2/r)) * log(k);
% m_range = round(core * [1]);
% for m=m_range
%    % computing p_tilde
%    Xtilde = uniform_sampling(X_skewed, m); % subsample X
%    h_tilde = efficient_H_diagonal(Xtilde);
%    p_tilde = h_tilde / k;
%    % test uniformity
%    pd = makedist('Multinomial', 'probabilities', p_tilde);
%    idx_sample = random(pd, m, 1);
%    %figure
%    %cdfplot(idx_sample)
%    test_pd = makedist('Uniform', 'Lower', 1, 'Upper', m);
%    [hypothesis, p_value, test_stat] = kstest(idx_sample, 'CDF', test_pd);
%    fprintf('\n --- m = %d', m)   
%    fprintf('\n - Uniformity test p-value: %.4f', p_value)
%    sum_pvalue_tilde = sum_pvalue_tilde + p_value;
% end
% fprintf('\n')
% end

%% C(iii) revised

rng(20190707)
% data matrices
n = 10^6;
k = 10;
pd_skewed = makedist('Gamma', 'a', 5, 'b', 14);
pd_normal = makedist('Normal');
X_normal = random(pd_normal, n, k);  % N(0,1) data
X_skewed = random(pd_skewed, n, k);  % skewed data

% resampling size
r = 4;
m_lower = round((n ^ (2/r)) * (k^(1+2/r)) * log(k));

% Simulation
reps = 500;
for m = m_lower * [1, 2]
    fprintf('\n ======== m = %d', m)
    [normdiff_normal_leverage, normdiff_normal_uniform, ...
     normdiff_skewed_leverage, normdiff_skewed_uniform] = ...
        q2c_p2_simulation(X_normal, X_skewed, m, reps);
    delta_normal = normdiff_normal_uniform - normdiff_normal_leverage;
    delta_skewed = normdiff_skewed_uniform - normdiff_skewed_leverage;
    fprintf('\n ---- Normal: ')
    fprintf('\n Mean Difference: %.10f', mean(delta_normal))
    fprintf('\n Variance in the difference: %10f', var(delta_normal))
    fprintf('\n Mean Squared Difference(1e-6): %.10f', mean(delta_normal.^2)/10^6)
    fprintf('\n Mean Absolute Difference(1e-6): %.10f', mean(abs(delta_normal))/10^6)
    fprintf('\n ---- Skewed: ')
    fprintf('\n Mean Difference: %.10f', mean(delta_skewed))
    fprintf('\n Variance in the difference: %10f', var(delta_skewed))
    fprintf('\n Mean Squared Difference(1e-6): %.10f', mean(delta_skewed.^2)/10^6)
    fprintf('\n Mean Absolute Difference(1e-6): %.10f', mean(abs(delta_skewed))/10^6)
    fprintf('\n ---- Testing H_0: mean_distortion_skewed>mean_distortion_normal')
    fprintf('\n z-stat: %.10f', sqrt(reps) *...
        (mean(delta_skewed)-mean(delta_normal))/...
        (std(delta_normal)+std(delta_skewed)))
    fprintf('\n')
end