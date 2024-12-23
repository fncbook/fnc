  function [x,u] = bvp(phi, a, b, ga, gb, init)
%BVP      Solve a boundary-value problem by finite differences
%         with either Dirichlet or Neumann BCs.
% Input:
%   phi      defines u'' = phi(x,u,u') (function)
%   a, b     endpoints of the domain (scalars)
%   ga       residual boundary function of u(a), u'(a) 
%   gb       residual boundary function of u(b), u'(b) 
%   init     initial guess for the solution (length n+1 vector)
% Output:
%   x        nodes in x (vector, length n+1)
%   u        values of u(x)  (vector, length n+1)
%   res      function for computing the residual

n = length(init) - 1;
[x, Dx, Dxx] = diffmat2(n, [a, b]);
h = x(2) - x(1);

u = levenberg(@residual, init);
u = u(:, end);

    function f = residual(u)
        % Computes the difference between u'' and phi(x,u,u') at the
        % interior nodes and appends the error at the boundaries. 
        du_dx = Dx*u;                   % discrete u'
        d2u_dx2 = Dxx*u;                % discrete u''
        f = d2u_dx2 - phi(x,u,du_dx);
        
        % Replace first and last values by boundary conditions.
        f(1) =   ga(  u(1),   du_dx(1)) / h;
        f(end) = gb(u(end), du_dx(end)) / h;
    end

end