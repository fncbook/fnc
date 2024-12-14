---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---
# Chapter 6

## Functions

(function-euler-matlab)=
``````{dropdown} Euler's method for an initial-value problem
```{literalinclude} ../matlab/fnc/euler.m
:language: matlab
:linenos: true
```
::::{admonition} About the code
:class: dropdown
The structure created by `ODEFunction` contains the data used to define it, and it is accessed in {numref}`Function {number} <function-euler>` by dot notation, such as `ivp.f`.
::::
``````

(function-ie2-matlab)=
``````{dropdown} Improved Euler method for an IVP
```{literalinclude} ../matlab/fnc/ie2.m
:language: matlab
:linenos: true
```
``````

(function-rk4-matlab)=
``````{dropdown} Fourth-order Runge-Kutta for an IVP
```{literalinclude} ../matlab/fnc/rk4.m
:language: matlab
:linenos: true
```
``````

(function-rk23-matlab)=
``````{dropdown} Adaptive IVP solver based on embedded RK formulas
```{literalinclude} ../matlab/fnc/rk23.m
:language: matlab
:linenos: true
```
::::{admonition} About the code
:class: dropdown
The check `t[i]+h==t[i]`on line 19 is to detect when $h$ has become so small that it no longer changes the floating-point value of $t_i$. This may be a sign that the underlying exact solution has a singularity near $t=t_i$, but in any case, the solver must halt by using a `break` statement to exit the loop.

On line 30, we use a combination of absolute and relative tolerances to judge the acceptability of a solution value, as in {eq}`absreltolerance`. In lines 41--43 we underestimate the step factor $q$ a bit and prevent a huge increase in the step size, since a rejected step is expensive, and then we make sure that our final step doesn't take us past the end of the domain.

Finally, line 37 exploits a subtle property of the BS23 formula called *first same as last* (FSAL).
While {eq}`bs23` calls for four stages to find the paired second- and third-order estimates, the final stage computed in stepping from $t_i$ to $t_{i+1}$ is identical to the first stage needed to step from $t_{i+1}$ to $t_{i+2}$. By repurposing `s₄` as `s₁` for the next pass, one of the stage evaluations comes for free, and only three evaluations of $f$ are needed per successful step.
::::

``````

(function-ab4-matlab)=
``````{dropdown} 4th-order Adams–Bashforth formula for an IVP
```{literalinclude} ../matlab/fnc/ab4.m
:language: matlab
:linenos: true
```
::::{admonition} About the code
:class: dropdown
Line 15 sets `σ` to be the coefficients of the generating polynomial $\sigma(z)$ of AB4. Lines 19--21 set up the IVP over the time interval $a \le t \le a+3 h$, call `rk4` to solve it using the step size $h$, and use the result to fill the first four values of the solution. Then line 24 computes the vector $[f_2,f_1,f_0]$.

Line 28 computes $f_i$, based on the most recent solution value and time. That goes into the first spot of `f`, followed by the three values that were previously most recent. These are the four values that appear in {eq}`ab4`. Each particular $f_i$ value starts at the front of `f`, moves through each position in the vector over three iterations, and then is forgotten.
::::
``````

(function-am2-matlab)=
``````{dropdown} 2nd-order Adams–Moulton (trapezoid) formula for an IVP
```{literalinclude} ../matlab/fnc/am2.m
:language: matlab
:linenos: true
```
::::{admonition} About the code
:class: dropdown
Lines 22-23 define the function $\mathbf{g}$ and call `levenberg` to find the new solution value, using an Euler half-step as its starting value. A robust code would have to intercept the case where `levenberg` fails to converge, but we have ignored this issue for the sake of brevity.
::::
``````

## Examples

```{code-cell}
:tags: [remove-cell]
addpath /Users/driscoll/Documents/GitHub/fnc/matlab/
addpath /Users/driscoll/Documents/GitHub/fnc/matlab/fnc
format short
set(0, 'defaultaxesfontsize', 12)
set(0, 'defaultlinelinewidth', 1.5, 'defaultscattermarkerfacecolor', 'k')
set(0, 'defaultFunctionLinelinewidth', 1.5)
```
### Section 6.1

(demo-basics-first-matlab)=
``````{dropdown} Solving an IVP
Let's use it to define and solve an initial-value problem for $u'=\sin[(u+t)^2]$ over $t \in [0,4]$, such that $u(0)=-1$.

To create an initial-value problem for $u(t)$, you must create an `ode` with a function that computes $u'$ and an initial condition for $u$. Then you create a solution by calling `solve` with a time interval. 

```{index} ! MATLAB; ode, ! MATLAB; solve
```

```{code-cell}
ivp = ode;
ivp.ODEFcn = @(t,u) sin((t + u)^2);
ivp.InitialTime = 0;
ivp.InitialValue = -1;
sol = solve(ivp, 0, 4);
```

The resulting solution object has fields `Time` and `Solution` that contain the approximate values of the solution at automatically chosen times in the interval you provided.

```{code-cell}
clf
plot(sol.Time, sol.Solution, '-o')
xlabel("t"), ylabel("u(t)")
title("Solution of an IVP")
```

You might want to know the solution at particular times other than the ones selected by the solver. That requires an interpolation, which is done by `solutionFcn`.

```{code-cell}
u = solutionFcn(ivp, 0, 10);
u(0:5)
```
``````

(demo-basics-sing-matlab)=
``````{dropdown} Finite-time singularity

The equation $u'=(u+t)^2$ gives us some trouble.

```{code-cell}
ivp = ode;
ivp.ODEFcn = @(t,u) (t + u)^2;
ivp.InitialTime = 0;
ivp.InitialValue = 1;
sol = solve(ivp, 0, 1);
```

The warning message we received can mean that there is a bug in the formulation of the problem. But if everything has been done correctly, it suggests that the solution may not exist past the indicated time. This is a possibility in nonlinear ODEs.

```{code-cell}
clf
semilogy(sol.Time, sol.Solution,)
xlabel("t"), ylabel("u(t)")
title("Finite-time blowup")
```
``````

(demo-basics-cond-matlab)=
``````{dropdown} Conditioning of an IVP
Consider the ODEs $u'=u$ and $u'=-u$. In each case we compute $\partial f/\partial u = \pm 1$, so the condition number bound from {numref}`Theorem %s <theorem-depIC>` is $e^{b-a}$ in both problems. However, they behave quite differently. In the case of exponential growth, $u'=u$, the bound is the actual condition number.

```{code-cell}
clf
for u0 = [0.7, 1, 1.3]    % initial values
    fplot(@(t) exp(t) * u0, [0, 3]), hold on
end
xlabel('t'), ylabel('u(t)')   % ignore this line
title('Exponential divergence of solutions')   % ignore this line
```

But with $u'=-u$, solutions actually get closer together with time.

```{code-cell}
clf
for u0 = [0.7, 1, 1.3]    % initial values
    fplot(@(t) exp(-t) * u0, [0, 3]), hold on
end
xlabel('t'), ylabel('u(t)')   % ignore this line
title('Exponential convergence of solutions')   % ignore this line
```

In this case the actual condition number is one, because the initial difference between solutions is the largest over all time. Hence the exponentially growing bound $e^{b-a}$ is a gross overestimate.
``````

### Section 6.2
(demo-euler-converge-matlab)=
``````{dropdown} Convergence of Euler's method
We consider the IVP $u'=\sin[(u+t)^2]$ over $0\le t \le 4$, with $u(0)=-1$. We need to define the function for the right-hand side of the ODE, the interval for the independent variable, and the initial value.

```{code-cell}
ivp = ode;
ivp.ODEFcn = @(t, u) sin((t + u)^2);
ivp.InitialTime = 0;
ivp.InitialValue = -1;
```

Here is the call to {numref}`Function {number} <function-euler>`.

```{code-cell}
[t, u] = eulerivp(ivp, b, 20);
clf, plot(t, u, '.-')
xlabel('t'), ylabel('u(t)')
title('Solution by Euler''s method')  % ignore this line
```

We could define a different interpolant to get a smoother picture above, but the derivation of Euler's method assumed a piecewise linear interpolant. We can instead request more steps to make the interpolant look smoother.

```{code-cell}
[t, u] = eulerivp(ivp, b, 50);
hold on, plot(t, u, '.-')
legend('20 steps', '50 steps')
```

Increasing $n$ changed the solution noticeably. Since we know that interpolants and finite differences become more accurate as $h\to 0$, we should anticipate the same behavior from Euler's method. We don't have an exact solution to compare to, so we will use a built-in solver to construct an accurate reference solution.

```{code-cell}
ivp.AbsoluteTolerance = 1e-13;
ivp.RelativeTolerance = 1e-13;
u_exact = solutionFcn(ivp, a, b);
```

Now we can perform a convergence study.

```{code-cell}
n = round(5 * 10.^(0:0.5:3));
err = [];
for k = 1:length(n)
    [t, u] = eulerivp(ivp, b, n(k));
    err(k) = norm(u_exact(t) - u, Inf);
end
table(n', err', VariableNames=["n", "inf-norm error"])
```

The error is approximately cut by a factor of 10 for each increase in $n$ by the same factor. A log-log plot also confirms first-order convergence. Keep in mind that since $h=(b-a)/n$, it follows that $O(h)=O(n^{-1})$.

```{code-cell}
clf
loglog(n, err, 'o-')
hold on, loglog(n, 0.5 * err(end) * (n / n(end)).^(-1), '--')
xlabel('n'), ylabel('inf-norm error')  % ignore this line
title('Convergence of Euler''s method')  % ignore this line
legend('error', '$O(n^{-1})$', 'location', 'southwest')  % ignore this line
```
``````

### Section 6.3
(demo-systems-predator-matlab)=
``````{dropdown} Predator-prey model
We encode the predator–prey equations via a function, defined here externally.

```{literalinclude} f63_predprey.m
:language: matlab
```

In this case, the ODE function accepts the required inputs `t` and `u` as well as a vector of parameters whose values don't change throughout a single instance of the problem. Now we can specify and solve the IVP like before, also giving a value for the parameter vector.

```{code-cell}
u0 = [1; 0.01];    # column vector
p = [0.1, 0.25];
ivp = ode(
    ODEFcn=@f63_predprey,...
    InitialTime=0,...
    InitialValue=u0,...
    Parameters=p,...
    );
sol = solve(ivp, 60);
size(sol.Solution)
```

Each column of the `Solution` field is the solution vector $\mathbf{u}$ at a particular time; each row is a component of $\mathbf{u}$ over all time.

```{code-cell}
clf
plot(sol.Time, sol.Solution)
xlabel('t'), ylabel('u(t)'), title('Predator-prey solution')  % ignore this line
legend('prey', 'predator')  % ignore this line
```

We can also use {numref}`Function {number} <function-euler>` to find the solution.

```{code-cell}
[t, u] = eulersys(ivp, 60, 1200);
```

```{code-cell}
hold on
plot(t, u, '.')
```

Notice above that the accuracy of the Euler solution deteriorates rapidly.

When there are just two components, it's common to plot the solution in the _phase plane_, i.e., with $u_1$ and $u_2$ along the axes and time as a parameterization of the curve.

```{code-cell}
clf
plot(u(1, :), u(2, :))
title("Predator-prey in the phase plane")
xlabel("y"), ylabel=("z")
```

From this plot we can deduce that the solution approaches a periodic one, which in the phase plane is represented by a closed loop.
``````

(demo-systems-coupledpendula-matlab)=
``````{dropdown} Coupled pendulums
::::{grid} 1 1 2 2
Let's implement the coupled pendulums from {numref}`Example {number} <example-systems-coupledpendula>`. The pendulums will be pulled in opposite directions and then released together from rest.
:::{card}
The `similar` function creates an array of the same size and type as a given value, without initializing the contents.
:::
::::

```{literalinclude} f63_pendulums.m
:language: matlab
```


```{code-cell}
function couple(u, p, t)
    γ, L, k = p
    g = 9.8
    udot = similar(u)
    udot[1:2] .= u[3:4]
    udot[3] = -γ * u[3] - (g / L) * sin(u[1]) + k * (u[2] - u[1])
    udot[4] = -γ * u[4] - (g / L) * sin(u[2]) + k * (u[1] - u[2])
    return udot
end

u₀ = [1.25, -0.5, 0, 0]
tspan = (0.0, 50.0);
```

::::{grid} 1 1 2 2
First we check the behavior of the system when the pendulums are uncoupled, i.e., when $k=0$.
:::{card}
Here `idxs` is used to plot two components as functions of time.
:::
::::

```{code-cell}
γ, L, k = 0, 0.5, 0
ivp = ODEProblem(couple, u₀, tspan, [γ, L, k])
sol = solve(ivp, Tsit5())
plot(sol, idxs=[1, 2], label=[L"\theta_1" L"\theta_2"],
    xlims=[20, 50], title="Uncoupled pendulums")
```

You can see that the pendulums swing independently. Because the model is nonlinear and the initial angles are not small, they have slightly different periods of oscillation, and they go in and out of phase.

With coupling activated, a different behavior is seen.

```{code-cell}
k = 1
ivp = ODEProblem(couple, u₀, tspan, [γ, L, k])
sol = solve(ivp, Tsit5())
plot(sol, idxs=[1, 2], label=[L"\theta_1" L"\theta_2"],
    xlims=[20, 50], title="Coupled pendulums")
```

The coupling makes the pendulums swap energy back and forth.
``````

(demo-rk-converge-matlab)=
``````{dropdown} Convergence of Runge–Kutta methods
We solve the IVP $u'=\sin[(u+t)^2]$ over $0\le t \le 4$, with $u(0)=-1$.

```{code-cell}
f(u, p, t) = sin((t + u)^2)
tspan = (0.0, 4.0)
u₀ = -1.0

ivp = ODEProblem(f, u₀, tspan)
```

We use a `DifferentialEquations` solver to construct an accurate approximation to the exact solution.

```{code-cell}
u_ref = solve(ivp, Tsit5(), reltol=1e-14, abstol=1e-14);
```

Now we perform a convergence study of our two Runge–Kutta implementations.

```{code-cell}
n = [round(Int, 2 * 10^k) for k in 0:0.5:3]
err_IE2, err_RK4 = [], []
for n in n
    t, u = FNC.ie2(ivp, n)
    push!(err_IE2, maximum(@.abs(u_ref(t) - u)))
    t, u = FNC.rk4(ivp, n)
    push!(err_RK4, maximum(@.abs(u_ref(t) - u)))
end

pretty_table([n err_IE2 err_RK4], header=["n", "IE2 error", "RK4 error"])
```

The amount of computational work at each time step is assumed to be proportional to the number of stages. Let's compare on an apples-to-apples basis by using the number of $f$-evaluations on the horizontal axis.

```{code-cell}
plot([2n 4n], [err_IE2 err_RK4], m=3, label=["IE2" "RK4"],
    xaxis=(:log10, "f-evaluations"), yaxis=(:log10, "inf-norm error"),
    title="Convergence of RK methods", leg=:bottomleft)

plot!(2n, 1e-5 * (n / n[end]) .^ (-2), l=:dash, label=L"O(n^{-2})")
plot!(4n, 1e-10 * (n / n[end]) .^ (-4), l=:dash, label=L"O(n^{-4})")
```

The fourth-order variant is more efficient in this problem over a wide range of accuracy.
``````

(demo-adapt-basic-matlab)=
``````{dropdown} Adaptive step size
Let's run adaptive RK on  $u'=e^{t-u\sin u}$.

```{code-cell}
f(u, p, t) = exp(t - u * sin(u))
ivp = ODEProblem(f, 0, (0.0, 5.0))
t, u = FNC.rk23(ivp, 1e-5)

plot(t, u, m=2,
    xlabel=L"t", ylabel=L"u(t)", title="Adaptive IVP solution")
```

The solution makes a very abrupt change near $t=2.4$. The resulting time steps vary over three orders of magnitude.

```{code-cell}
Δt = diff(t)
plot(t[1:end-1], Δt, title="Adaptive step sizes",
    xaxis=(L"t", (0, 5)), yaxis=(:log10, "step size"))
```

If we had to run with a uniform step size to get this accuracy, it would be

```{code-cell}
println("minimum step size = $(minimum(Δt))")
```

On the other hand, the average step size that was actually taken was

```{code-cell}
println("average step size = $(sum(Δt)/(length(t)-1))")
```

We took fewer steps by a factor of almost 1000! Even accounting for the extra stage per step and the occasional rejected step, the savings are clear.

``````

(demo-adpat-sing-matlab)=
``````{dropdown} Adaptive step size near a singularity
In {numref}`Demo %s <demo-basics-sing>` we saw an IVP that appears to blow up in a finite amount of time. Because the solution increases so rapidly as it approaches the blowup, adaptive stepping is required even to get close.

```{code-cell}
f(u, p, t) = (t + u)^2
ivp = ODEProblem(f, 1, (0.0, 1.0))
t, u = FNC.rk23(ivp, 1e-5);
```

In fact, the failure of the adaptivity gives a decent idea of when the singularity occurs.

```{code-cell}
plot(t, u, legend=:none,
    xlabel=L"t", yaxis=(:log10, L"u(t)"), title="Finite-time blowup")

tf = t[end]
vline!([tf], l=:dash)
annotate!(tf, 1e5, latexstring(@sprintf("t = %.6f ", tf)), :right)
```
``````

(demo-implicit-ab4-matlab)=
``````{dropdown} Convergence of Adams–Bashforth
We study the convergence of AB4 using the IVP $u'=\sin[(u+t)^2]$ over $0\le t \le 4$, with $u(0)=-1$. As usual, `solve` is called to give an accurate reference solution.

```{code-cell}
ivp = ODEProblem((u, p, t) -> sin((t + u)^2), -1.0, (0.0, 4.0))
u_ref = solve(ivp, Tsit5(), reltol=1e-14, abstol=1e-14);
```

Now we perform a convergence study of the AB4 code.

```{code-cell}
n = @. [round(Int, 4 * 10^k) for k in 0:0.5:3]
err = []
for n in n
    t, u = FNC.ab4(ivp, n)
    push!(err, norm(u_ref.(t) - u, Inf))
end

pretty_table([n err], header=["n", "inf-norm error"])
```

The method should converge as $O(h^4)$, so a log-log scale is appropriate for the errors.

```{code-cell}
plot(n, err, m=3, label="AB4",
    xaxis=(:log10, L"n"), yaxis=(:log10, "inf-norm error"),
    title="Convergence of AB4", leg=:bottomleft)

plot!(n, (n / n[1]) .^ (-4), l=:dash, label=L"O(n^{-4})")
```
``````

(demo-implicit-stiff-matlab)=
``````{dropdown} Stiffness
The following simple ODE uncovers a surprise.

```{code-cell}
ivp = ODEProblem((u, p, t) -> u^2 - u^3, 0.005, (0, 400.0))
```

We will solve the problem first with the implicit AM2 method using $n=200$ steps.

```{code-cell}
tI, uI = FNC.am2(ivp, 200)

plot(tI, uI, label="AM2",
    xlabel=L"t", ylabel=L"u(t)", leg=:bottomright)
```

Now we repeat the process using the explicit AB4 method.

```{code-cell}
tE, uE = FNC.ab4(ivp, 200)

scatter!(tE, uE, m=3, label="AB4", ylim=[-4, 2])
```

Once the solution starts to take off, the AB4 result goes catastrophically wrong.

```{code-cell}
uE[105:111]
```

We hope that AB4 will converge in the limit $h\to 0$, so let's try using more steps.

```{code-cell}
plt = scatter(tI, uI, label="AM2, n=200", m=3,
    xlabel=L"t", ylabel=L"u(t)", leg=:bottomright)

for n in [1000, 1600]
    tE, uE = FNC.ab4(ivp, n)
    plot!(tE, uE, label="AM4, n=$n")
end
plt
```

So AB4, which is supposed to be _more_ accurate than AM2, actually needs something like 8 times as many steps to get a reasonable-looking answer!
``````

(demo-zs-LIAF-matlab)=
``````{dropdown} Instability
We'll measure the error at the time $t=1$.

```{code-cell}
du_dt(u, t) = u
û = exp
a, b = 0.0, 1.0;
n = [5, 10, 20, 40, 60]
err = []
t, u = [], []
for n in n
    h = (b - a) / n
    t = [a + i * h for i in 0:n]
    u = [1; û(h); zeros(n - 1)]
    f_val = [du_dt(u[1], t[1]); zeros(n)]
    for i in 2:n
        f_val[i] = du_dt(u[i], t[i])
        u[i+1] = -4 * u[i] + 5 * u[i-1] + h * (4 * f_val[i] + 2 * f_val[i-1])
    end
    push!(err, abs(û(b) - u[end]))
end

pretty_table([n (b - a) ./ n err], header=["n", "h", "error"])
```

The error starts out promisingly, but things explode from there. A graph of the last numerical attempt yields a clue.

```{code-cell}
plot(t, abs.(u), m=3, label="",
    xlabel=L"t", yaxis=(:log10, L"|u(t)|"), title="LIAF solution")
```

It's clear that the solution is growing exponentially in time.
``````
