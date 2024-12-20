function [I, x] = intinf(f, tol)
% INTINF   Adaptive doubly exponential integration over (-inf,inf).
% Input:
%   f   integrand (function)
%   tol error tolerance (positive scalar)
% Output:
%   I   approximation to intergal(f) over (-inf,inf)
%   x   evaluation nodes (vector) 

xi = @(t) sinh(sinh(t));
dxi_dt = @(t) cosh(t) .* cosh(sinh(t));
g = @(t) f(xi(t)) .* dxi_dt(t);

% Find where to truncate the integration interval.
M = 3;
while (abs(g(-M)) > tol/100) || (abs(g(M)) > tol/100)
    M = M + 0.5;
    if isinf(x(M)) 
        warning("Function may not decay fast enough.")
        M = M - 0.5;
        break
    end
end

[I, t] = intadapt(g, -M, M, tol);
x = xi(t);
end