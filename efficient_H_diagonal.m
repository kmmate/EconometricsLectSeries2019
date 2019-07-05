function h = efficient_H_diagonal(X)
%% Efficient computation of the diagonal elements of H, according to
% Algorithm 1 in the latex file
%%
n = size(X, 1);
B = inv(X' * X);
h = zeros(n,1);
for i=1:n
   h(i) = X(i,:) * B * X(i,:)'; 
end
end