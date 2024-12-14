function [t, u] = eulerivp(ivp, a, b, n)
% EULERIVP   Euler's method for a scalar initial-value problem.
% Input:
%   ivp     structure defining the IVP
%   a, b    endpoints of time interval (scalars)
%   n       number of time steps (integer)
% Output:
%   t       selected nodes (vector, length n+1)
%   u       solution values (array, m by (n+1))

du_dt = ivp.ODEFcn;
u0 = ivp.InitialValue;
p = ivp.Parameters;

% Define time discretization.
h = (b - a) / n;
t = a + (0:n) * h;

% Initialize solution array.
u = zeros(length(u0), n+1);
u(:, 1) = u0;

% Time stepping.
for i = 1:n
  u(:, i+1) = u(:, i) + h * du_dt(t(i), u(i), p);
end