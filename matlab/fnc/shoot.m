function [x, u, du_dx] = shoot(phi, a, b, ga, gb, init, tol)
%SHOOT   Shooting method for a two-point boundary-value problem.
% Input:
%   phi      defines u'' = phi(x, u, u') (function)
%   a, b     endpoints of the domain (scalars)
%   ga       residual boundary function of u(a), u'(a) 
%   gb       residual boundary function of u(b), u'(b) 
%   init     initial guess for u(a) and u'(a) (column vector)
%   tol      error tolerance (scalar)
% Output:
%   x        nodes in x (length n+1)
%   u        values of u(x)  (length n+1)
%   du_dx    values of u'(x) (length n+1)

% To be solved by the IVP solver
function f = shootivp(x, y)
  f = [ y(2); phi(x, y(1), y(2)) ];
end

ivp = ode(ODEFcn=@shootivp);
ivp.InitialTime = a;
ivp.AbsoluteTolerance = tol / 10;
ivp.RelativeTolerance = tol / 10;

% To be solved by levenberg
function residual = objective(s)
  ivp.InitialValue = s;
  sol = solve(ivp, a, b);
  x = sol.Time;  y = sol.Solution;
  residual = [ga(y(1, 1), y(2, 1)); gb(y(1, end), y(2, end))];
end

y = [];    % shared variable
s = levenberg(@objective, init, tol);

% Don't need to solve the IVP again. It was done within the
% objective function already.
u = y(1, :);            % solution     
du_dx = y(2, :);         % derivative

end