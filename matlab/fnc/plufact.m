function [L, U, p] = plufact(A)
% PLUFACT   Pivoted LU factorization 
% Input:
%   A    square matrix
% Output:
%   L    unit lower triangular 
%   U    upper triangular
%   p    row permutation vector such that A(p,:) = L*U
n = size(A, 1);
L = zeros(n, n);
U = zeros(n, n);
p = zeros(1, n);
A_k = A;

% Reduction by outer products
for k = 1:n
    [~, p(k)] = max(abs(A_k(:, k)));
    U(k, :) = A_k(p(k), :);
    L(:, k) = A_k(:, k) / U(k, k);
    if k < n
        A_k = A_k - L(:, k) * U(k, :);
    end
end
L = tril(L(p, :));
U = triu(U);
end
