function [t, u] = rk23(ivp, a, b, tol)
% RK23   Adaptive IVP solver based on embedded RK formulas.
% Input:
%   ivp     structure defining the IVP
%   a, b    endpoints of time interval (scalars)
%   tol     global error target (positive scalar)
% Output:
%   t       selected nodes (vector, length n+1)
%   u       solution values (array, m by (n+1))

du_dt = ivp.ODEFcn;
u0 = ivp.InitialValue;
p = ivp.Parameters;

% Initialize for the first time step.
t = a;
u(:, 1) = u0(:);  i = 1;
h = 0.5 * tol^(1/3);
s1 = du_dt(t(1), u(:, 1));

% Time stepping.
while t(i) < b
    % Detect underflow of the step size.
    if t(i) + h == t(i)
        warning('Stepsize too small near t=%.6g.',t(i))
        break  % quit time stepping loop
    end
    
    % New RK stages.
    s2 = du_dt(t(i) + h/2,   u(:, i) + (h/2)   * s1, p);
    s3 = du_dt(t(i) + 3*h/4, u(:, i) + (3*h/4) * s2, p);
    unew2 = u(:, i) + h * (2*s1 + 3*s2 + 4*s3) / 9;    % 2rd order solution
    s4 = du_dt(t(i) + h, unew2, p );
    err = h * (-5*s1/72 + s2/12 + s3/9 - s4/8);        % 2nd/3rd order difference
    E = norm(err, Inf);                                % error estimate
    maxerr = tol * (1 + norm(u(:, i), Inf));           % relative/absolute blend
    
    % Accept the proposed step? 
    if E < maxerr     % yes 
        t(i+1) = t(i) + h;
        u(:, i+1) = unew2;
        i = i+1;
        s1 = s4;      % use FSAL property
    end
    
    % Adjust step size. 
    q = 0.8 * (maxerr/E)^(1/3);       % conservative optimal step factor
    q = min(q, 4);                    % limit stepsize growth
    h = min(q*h, b - t(i));           % don't step past the end
end