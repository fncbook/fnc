function [t, u] = ie2(ivp, a, b, n)
% IE2    Improved Euler method for an IVP.
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
    uhalf = u(:, i) + h/2 * du_dt(t(i), u(:, i), p);
    u(:, i+1) = u(:, i) + h * du_dt(t(i) + h/2, uhalf, p);
end