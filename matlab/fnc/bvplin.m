function [x, u] = bvplin(p, q, r, a, b, ua, ub,n)
% BVPLIN   Solve a linear boundary-value problem.
% Input:
%   p, q, r  u'' + pu' + qu = r (functions)
%   a, b     endpoints of the domain (scalars)
%   ua       value of u(a)
%   ub       value of u(b)
%   n        number of subintervals (integer)
% Output:
%   x       collocation nodes (vector, length n+1)
%   u       solution at nodes (vector, length n+1)

[x, Dx, Dxx] = diffmat2(n, [a, b]);

P = diag(p(x));
Q = diag(q(x));
L = Dxx + P * Dx + Q;    % ODE expressed at the nodes
r = r(x);    

% Replace first and last rows using boundary conditions. 
I = speye(n+1); 
A = [ I(:, 1)'; L(2:n, :); I(:, n+1)' ];

f = [ ua; r(2:n); ub ];

% Solve the system.
u = A \ f;