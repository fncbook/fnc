function [x, u] = parabolic(phi, xspan, m, ga, gb, tspan, init)
% PARABOLIC   Solve parabolic PDE by the method of lines.
% Input:
%   phi      defines ∂u/∂t = phi(t, x, u, ∂u/∂x, ∂^2u/∂x^2) 
%   xspan    spatial domain
%   m        number of spatial nodes
%   ga, gb   boundary conditions as functions of u and ∂u/∂x
%   tspan    time interval
%   init     initial condition as a function of x
% Output       
%   x        spatial nodes (vector)
%   u        function for the solution u(t) at nodes

    [x, Dx, Dxx] = diffcheb(m, xspan);
    int = 2:m;    % indexes of interior nodes

    function u = extend(v)
        function residual = objective(ubc)
            ua = ubc(1);  ub = ubc(2);
            ux = Dx * [ua; v; ub];
            residual = [ga(ua, ux(1)); gb(ub, ux(end))];
        end
        ubc = levenberg(@objective, [0, 0]);
        ubc = ubc(:, end);
        u = [ubc(1); v; ubc(2)];
    end

    function f = mol_ode(t, v, p)
        u = extend(v);
        ux = Dx * u;
        uxx = Dxx * u;
        f = phi(t, x(int), u(int), ux(int), uxx(int));
    end

    ivp = ode(ODEFcn=@mol_ode);
    ivp.InitialTime = tspan(1);
    ivp.InitialValue = init(x(int));
    ivp.Solver = "stiff";
    sol = solutionFcn(ivp, tspan(1), tspan(2));
    u = @(t) extend(sol(t));
end