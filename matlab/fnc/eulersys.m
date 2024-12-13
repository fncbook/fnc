function [t,u] = eulersys(ivp, a, b, n)
% EULERSYS   Euler's method for a first-order IVP system.
% Input:
%   du_dt   defines f in u'(t)=f(t,u) (function)
%   tspan   endpoints of time interval (2-vector)
%   u0      initial value (vector, length m)
%   n       number of time steps (integer)
% Output:
%   t       selected nodes  (vector, length n+1)
%   u       solution values   (array, (n+1)-by-m)

du_dt = ivp.ODEFcn;
u0 = ivp.InitialValue;
p = ivp.Parameters;
h = (b - a) / n;
t = a + (0:n) * h;
u = zeros(length(u0), n+1);
u(:, 1) = u0(:);
for i = 1:n
  u(:, i+1) = u(:, i) + h * du_dt(t(i), u(:, i), p);
end
