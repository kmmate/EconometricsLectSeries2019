function subX = uniform_sampling(X, m)
%% Sampling uniformly m rows of X
[n, k] = size(X);
%subX = zeros(m, k);
pd = makedist('Multinomial', 'probabilities', ones(n,1)/n);
rows = random(pd, m, 1);
subX = X(rows, :);
end