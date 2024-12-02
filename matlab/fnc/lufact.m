function [L, U] = lufact(A)
% LUFACT   LU factorization (demo only--not stable!).
% Input:
%   A    square matrix
% Output:
%   L,U  unit lower triangular and upper triangular such that LU=A

n = size(A, 1);     % detect the dimensions from the input
L = eye(n);         % ones on main diagonal, zeros elsewhere
U = zeros(n, n);
A_k = A;            % make a working copy 

% Reduction by outer products
for k = 1:n-1
    U(k, :) = A_k(k, :);
    L(:, k) = A_k(:,k) / U(k, k);
    A_k = A_k -  L(:, k) * U(k, :);
end
U(n, n) = A_k(n, n);
L = tril(L); U = triu(U);    % force exact triangularity