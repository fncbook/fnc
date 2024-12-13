function [t, u] = ab4(ivp, a, b, n)
%AB4     4th-order Adams-Bashforth formula for an IVP.
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

% Constants in the AB4 method.
k = 4;  
sigma = [55; -59; 37; -9] / 24;  

% Find starting values by RK4.
[ts, us] = rk4(ivp, a, a + (k-1)*h, k-1);
u = zeros(length(u0), n+1);
u(:, 1:k) = us(:, 1:k);

% Compute history of u' values, from oldest to newest.
f = zeros(length(u0), k);
for i = 1:k-1
  f(:, k-i) = du_dt(t(i), u(:, i), p);
end

% Time stepping.
for i = k:n
  f = [du_dt(t(i), u(:, i), p), f(:, 1:k-1)];   % new value of du/dt
  u(:, i+1) = u(:, i) + h * (f * sigma);        % advance one step
end