function [t, u] = am2(ivp, a, b, n)
% AM2    2nd-order Adams-Moulton (trapezoid) formula for an IVP.
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

u = zeros(length(u0), n+1);
u(:, 1) = u0(:);

% Time stepping.
for i = 1:n
  % Data that does not depend on the new value.
  known = u(:,i) + h/2 * du_dt(t(i), u(:, i), p);
  % Find a root for the new value. 
  unew = levenberg(@trapzero, known);
  u(:, i+1) = unew(:, end);
end

% This function defines the rootfinding problem at each step.
function F = trapzero(z)
    F = z - h/2 * du_dt(t(i+1), z, p) - known;
end

end  % main function