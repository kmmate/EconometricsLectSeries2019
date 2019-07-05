%& Tinbergen Econometrics Lecture Series 2019, Exam, Question 2 / a
%%
% given draws
h = [2 3 1 2 1 1 3 2 1];
g = [-1 -1 1 -1 1 -1 1 -1 1];
m = 3;  % nobs in sketched data
n = 9;  % nobs in original data
k = 4;  % no of variables
% construct Pi
P = zeros(m, n);
for i=1:n
   P(h(i),i) = g(i); 
end

% simulation
reps = 1000;
m_distance = @(A, B) norm(A-B,'fro');
distance_sum = 0.0;
for rep = 1:reps
   A = randn(n, k);
   A_tilde1 = P * A;
   A_tilde2 = countsketch(A, h, g, m);
   distance_sum = distance_sum + m_distance(A_tilde1, A_tilde2);
end
display(distance_sum / reps, ' Mean distance across reps ')

