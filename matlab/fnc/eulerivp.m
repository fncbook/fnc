function [t, u] = eulerivp(ivp, a, b, n)
% EULERIVP   Euler's method for a scalar initial-value problem.
% Input:
%   du_dt   defines f in u'(t)=f(t,u). (function)
%   tspan   endpoints of time interval (2-vector)
%   u0      initial value (scalar) 
%   n       number of time steps (integer)
% Output:
%   t       selected nodes  (vector, length n+1)
%   u       solution values   (vector, length n+1)

du_dt = ivp.ODEFcn;
u0 = ivp.InitialValue;
p = ivp.Parameters;
h = (b - a) / n;
t = a + (0:n) * h;
u = zeros(1, n+1);
u(1) = u0;
for i = 1:n
  u(i+1) = u(i) + h * du_dt(t(i), u(i), p);
end