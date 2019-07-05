function a_tilde = countsketch(A, h, g, m)
%% Streaming version of CountSketch algorithm
% Inputs
% ------
% - A: original data, matrix of size n-by-k
% - h: draws of h
% - g: draws of g
% - m: nobs in sketched data matrix

% Output
% ------
% - a_tilde: sketched data matrix of siye m-by-k

n = size(A, 1);
k = size(A, 2);
a_tilde = zeros(m, k);
for i = 1:n
    a_tilde(h(i), :) = a_tilde(h(i), :) + g(i) * A(i, :);
end
end