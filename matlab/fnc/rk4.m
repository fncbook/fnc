function [t, u] = rk4(ivp, a, b, n)
% RK4    Fourth-order Runge-Kutta for an IVP.
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
t = a + h * (0:n)';

% Initialize solution array. 
u = zeros(length(u0), n+1);
u(:, 1) = u0(:);

% Time stepping.
for i = 1:n
  k1 = h * du_dt( t(i),       u(:, i)       , p);
  k2 = h * du_dt( t(i) + h/2, u(:, i) + k1/2, p );
  k3 = h * du_dt( t(i) + h/2, u(:, i) + k2/2, p );
  k4 = h * du_dt( t(i) + h,   u(:, i) + k3  , p);
  u(:, i+1) = u(:, i) + (k1 + 2*(k2 + k3) + k4) / 6;
end