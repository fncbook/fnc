function [I, x] = intsing(f, tol)
% INTSING   Adaptively integrate a function with a singularity at the left endpoint.
% Input:
%   f   integrand  (function)
%   tol error tolerance (positive scalar)
% Output:
%   I   approximation to integral(f) over (0,1)
%   x   evaluation nodes (vector)

xi = @(t) 2 ./ (1 + exp( 2*sinh(t) ));
dxi_dt = @(t) cosh(t) ./ cosh( sinh(t) ).^2;
g = @(t) f(xi(t)) .* dxi_dt(t);

% Find where to truncate the integration interval.
M = 3;
while abs(g(M)) > tol/100
    M = M + 0.5;
    if iszero(x(M)) 
        warning("Function may grow too rapidly.")
        M = M - 0.5;
        break
    end
end

[I, t] = intadapt(g, 0, M, tol);
x = xi(t);
end