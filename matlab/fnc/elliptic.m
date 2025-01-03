function u = elliptic(phi, g, m, xspan, n, yspan)
%ELLIPTIC   Solve an elliptic PDE in 2d.
% Input:
%   phi          defines phi)(x,y,u,u_x,u_xx,u_y,u_yy) = 0 (function)
%   g            boundary condition (function)
%   m, xspan     size and interval of x discretization (integer, 2-vector)
%   n, yspan     size and interval of y discretization (integer, 2-vector)
% Output:
%   U            solution (n+1 by n+1)
%   X,Y          coordinate matrices (n+1 by n+1)

    % Discretize the domain.
    [x, Dx, Dxx] = diffcheb(m, xspan);
    [y, Dy, Dyy] = diffcheb(n, yspan);
    [mtx, X, Y, vec, unvec, is_boundary] = tensorgrid(x, y);

    % Identify boundary locations and evaluate the boundary condition.
    idx = vec(is_boundary);
    gb = g(X(idx), Y(idx));

    % Evaluate the PDE+BC residual.
    function r = residual(u)
        U = unvec(u);
        R = phi(X, Y, U, Dx * U, Dxx * U, U * Dy', U * Dyy');  % PDE
        R(idx) = u(idx) - gb;                                  % boundary
        r = vec(R);
    end

    % Solve the equation.
    u = levenberg(@residual, vec(zeros(size(X))));
    U = unvec(u(:, end));

    function u = evaluate(xi, eta)
        v = zeros(1, n+1);
        for j = 1:n+1
            v(j) = chebinterp(x, U(:, j), xi);
        end
        u = chebinterp(y, v, eta);
    end

    u = @evaluate;
end

function f = chebinterp(x, v, xi)
    n = length(x) - 1;
    w = (-1.0) .^ (0:n)';
    w([1, n+1]) = w([1, n+1]) / 2;

    terms = w ./ (xi - x(:));
    if any(isinf(terms))     % exactly at a node
        f = v(xi == x);
    else
        f = sum(v(:) .* terms) / sum(terms);
    end
end
