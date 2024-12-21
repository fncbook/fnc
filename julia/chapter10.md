---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
# Chapter 10

## Functions

(function-shoot-julia)=
``````{dropdown} Shooting method for a two-point boundary-value problem
```{literalinclude} ../julia/package/src/chapter10.jl
:filename: shoot.jl
:start-line: 0
:end-line: 29
:language: julia
:linenos: true
```

::::{admonition} About the code
:class: dropdown
Because `x` and `y` are assigned empty values in line 24, when the function `objective` runs it uses those values rather than new ones in local scope. Thus, line 19 updates them to hold the latest results of the IVP solver, saving the need to solve it again after `levenberg` has finished the rootfinding. 

The error tolerance in the IVP solver is kept smaller than in the rootfinder, to prevent the rootfinder from searching in a noisy landscape. Finally, note how line 28 uses destructuring of `eachrow(y)` to assign the columns of `y` to separate names.
::::

``````

(function-diffmats2-julia)=
``````{dropdown} Second-order differentiation matrices
```{literalinclude} ../julia/package/src/chapter10.jl
:filename: diffmats2.jl
:start-line: 31
:end-line: 62
:language: julia
:linenos: true
```
``````

(function-diffcheb-julia)=
``````{dropdown} Chebyshev differentiation matrices
```{literalinclude} ../julia/package/src/chapter10.jl
:filename: diffcheb.jl
:start-line: 64
:end-line: 92
:language: julia
:linenos: true
```
``````

(function-bvplin-julia)=
``````{dropdown} Solution of a linear boundary-value problem
```{literalinclude} ../julia/package/src/chapter10.jl
:filename: bvplin.jl
:start-line: 94
:end-line: 119
:language: julia
:linenos: true
```

::::{admonition} About the code
:class: dropdown
Note that there is no need to explicitly form the row-deletion matrix $\mathbf{E}$ from {eq}`rowdeletion`. Since it only appears as left-multiplying $\mathbf{L}$ or $\mathbf{r}$, we simply perform the row deletions as needed using indexing.
::::
``````

(function-bvp-julia)=
``````{dropdown} Solution of a nonlinear boundary-value problem
```{literalinclude} ../julia/package/src/chapter10.jl
:filename: bvp.jl
:start-line: 121
:end-line: 151
:language: julia
:linenos: true
```
:::{admonition} About the code
:class: dropdown
The nested function `residual` uses differentiation matrices computed externally to it, rather than computing them anew on each invocation. As in {numref}`Function {number} <function-bvplin>`, there is no need to form the row-deletion matrix $\mathbf{E}$ explicitly. In lines 23--24, we divide the values of $g_1$ and $g_2$ by a factor of $h$. This helps scale the residual components more uniformly and improves the robustness of convergence a bit.
:::
``````

(function-fem-julia)=
``````{dropdown} Piecewise linear finite elements for a linear BVP
```{literalinclude} ../julia/package/src/chapter10.jl
:filename: fem.jl
:start-line: 153
:end-line: 192
:language: julia
:linenos: true
```
``````

## Examples

```{code-cell}
:tags: [remove-output]
import Pkg; Pkg.activate("/Users/driscoll/Documents/GitHub/fnc")
using FundamentalsNumericalComputation
FNC.init_format()
```

### 10.1 @section-bvp-tpbvp

(demo-tpbvp-mems-julia)=
``````{dropdown} @demo-tpbvp-mems
```{index} ! Julia; in-place function
```

:::::{grid} 1 1 2 2
As a system, the MEMS problem from {numref}`Example {number} <example-tpbvp-mems>` uses $y_1=w$, $y_2=w'$ to obtain

:::{math}
:label: memssys-intro
\begin{split}
y_1' &= y_2, \\
y_2' &= \frac{\lambda}{y_1^2} - \frac{y_2}{r}.
\end{split}
:::

We will code an *in-place* form of this ODE, in which the first argument is used to return the computed values of $y_1'$ and $y_2'$.  
:::{card}
The in-place code here saves the computing time that would otherwise be needed to allocate memory for `f` repeatedly.
:::
:::::

```{code-cell}
function ode!(f, y, λ, r)
    f[1] = y[2]
    f[2] = λ / y[1]^2 - y[2] / r
    return nothing
end;
```

Notice that the `return` value is irrelevant with the in-place style. We use the same style for the boundary conditions $y_2(0)=0$, $y_1(2)=1$.

```{code-cell}
function bc!(g, y, λ, r)
    g[1] = y(0)[2]
    g[2] = y(1)[1] - 1
    return nothing
end;
```

In the `bc!` function, the `y` argument is just like an IVP solution from {numref}`section-ivp-basics`. Thus, `y(0)` is the value of the solution at $x=0$, and the second component of that value is what we wish to make zero. Similarly, `y(1)[1]` is the notation for $y_1(1)$, which is supposed to equal 1. 

The domain of the mathematical problem is $r\in [0,1]$. However, there is a division by $r$ in the ODE, so we want to avoid $r=0$ by truncating the domain a bit.

```{code-cell}
domain = (eps(), 1.0)
```

We need one last ingredient that is not part of the mathematical setup: an initial estimate for the solution. As we will see, this plays the same role as initialization in Newton's method for rootfinding. Here, we try a constant value for each component.

```{code-cell}
est = [1, 0]
```

Now we set up and solve a `BVProblem` with the parameter value $\lambda=0.6$.

```{code-cell}
bvp = BVProblem(ode!, bc!, est, domain, 0.6)
y = solve(bvp)
plot(
    y,
    label = [L"w" L"w'"],
    legend = :right,
    xlabel = L"r",
    ylabel = "solution",
    title = "Solution of MEMS problem for λ=0.6",
)
```

To visual accuracy, the boundary conditions have been enforced.
``````

### 10.2 @section-bvp-shooting
(demo-shooting-naive-julia)=
``````{dropdown} @demo-shooting-naive
::::{grid} 1 1 2 2

Let's first examine the shooting approach for the TPBVP from {numref}`Example {number} <example-tpbvp-mems>` with $\lambda=0.6$. 

:::{card}
The character `ϕ` is typed as `\phi`<kbd>Tab</kbd>. 
:::
::::

```{code-cell}
λ = 0.6
ϕ = (r, w, dwdr) -> λ / w^2 - dwdr / r;
```

We convert the ODE to a first-order system in order to apply a numerical method. We also have to truncate the domain to avoid division by zero.

```{code-cell}
f = (y, p, r) -> [y[2]; ϕ(r, y[1], y[2])]
a, b = eps(), 1.0;
```

The BVP specifies $w'(0)=y_2(0)=0$. We can try multiple values for the unknown $w(0)=y_1(0)$ and plot the solutions.

```{code-cell}
plt = plot(
    xaxis = (L"x"),
    yaxis = (L"w(x)"),
    title = "Different initial values",
    leg = :bottomright,
)

for w0 in 0.4:0.1:0.9
    IVP = ODEProblem(f, [w0, 0], (a, b))
    y = solve(IVP, Tsit5())
    plot!(y, idxs = [1], label = "w(0) = $w0")
end
plt
```

On the graph, it's the curve starting at $w(0)=0.8$ that comes closest to the required condition $w(1)=1$, but it's a bit too large.
``````

(demo-shooting-mems-julia)=
``````{dropdown} @demo-shooting-mems
We revisit {numref}`Demo {number} <demo-shooting-naive>` but let {numref}`Function {number} <function-shoot>` do the heavy lifting.

```{code-cell}
λ = 0.6
ϕ = (r, w, dwdr) -> λ / w^2 - dwdr / r;
a = eps();
b = 1;
```

We specify the given and unknown endpoint values.

```{code-cell}
g₁(w, dw) = dw     # w'=0 at left
g₂(w, dw) = w - 1    # w=1 at right

r, w, dw_dx = FNC.shoot(ϕ, (a, b), g₁, g₂, [0.8, 0])
plot(r, w, title = "Shooting solution", xaxis = (L"x"), yaxis = (L"w(x)"))
```

The value of $w$ at $r=1$, meant to be exactly one, was computed to be

```{code-cell}
@show w[end];
```

The accuracy is consistent with the error tolerance used for the IVP solution. The initial value $w(0)$ that gave this solution is

```{code-cell}
@show w[1];
```
``````

(demo-shooting-unstable-julia)=
``````{dropdown} @demo-shooting-unstable

```{code-cell}
plt = plot(
    xaxis = (L"x"),
    yaxis = ([-1.2, 0.5], L"u(x)"),
    title = "Shooting instability",
    leg = :topleft,
)
for λ in 6:4:18
    g₁(u, du) = u + 1
    g₂(u, du) = u
    ϕ = (x, u, du_dx) -> λ^2 * (u + 1)
    x, u = FNC.shoot(ϕ, (0.0, 1.0), g₁, g₂, [-1, 0])
    plot!(x, u, label = "λ=$λ")
end
plt
```

The numerical solutions evidently don't satisfy the right boundary condition as $\lambda$ increases, which makes them invalid. 

``````

### 10.3 @section-bvp-diffmats
(demo-diffmats-2nd-julia)=
``````{dropdown} @demo-diffmats-2nd
We test first-order and second-order differentiation matrices for the function $x + \exp(\sin 4x)$ over $[-1,1]$.

```{code-cell}
f = x -> x + exp(sin(4 * x));
```

For reference, here are the exact first and second derivatives.

```{code-cell}
dfdx = x -> 1 + 4 * exp(sin(4 * x)) * cos(4 * x);
d2fdx2 = x -> 4 * exp(sin(4 * x)) * (4 * cos(4 * x)^2 - 4 * sin(4 * x));
```

We discretize on equally spaced nodes and evaluate $f$ at the nodes.

```{code-cell}
t, Dₓ, Dₓₓ = FNC.diffmat2(18, [-1, 1])
y = f.(t);
```

Then the first two derivatives of $f$ each require one matrix-vector multiplication.

```{code-cell}
yₓ = Dₓ * y;
yₓₓ = Dₓₓ * y;
```

The results show poor accuracy for this small value of $n$.

```{code-cell}
plot(dfdx, -1, 1, layout = 2, xaxis = (L"x"), yaxis = (L"f'(x)"))
scatter!(t, yₓ, subplot = 1)
plot!(d2fdx2, -1, 1, subplot = 2, xaxis = (L"x"), yaxis = (L"f''(x)"))
scatter!(t, yₓₓ, subplot = 2)
```

An convergence experiment confirms the order of accuracy. Because we expect an algebraic convergence rate, we use a log-log plot of the errors.

```{code-cell}
n = @. round(Int, 2^(4:0.5:11))
err1 = zeros(size(n))
err2 = zeros(size(n))
for (k, n) in enumerate(n)
    t, Dₓ, Dₓₓ = FNC.diffmat2(n, [-1, 1])
    y = f.(t)
    err1[k] = norm(dfdx.(t) - Dₓ * y, Inf)
    err2[k] = norm(d2fdx2.(t) - Dₓₓ * y, Inf)
end
plot(n, [err1 err2], m = :o, label = [L"f'" L"f''"])
plot!(
    n,
    10 * 10 * n .^ (-2),
    l = (:dash, :gray),
    label = "2nd order",
    xaxis = (:log10, "n"),
    yaxis = (:log10, "max error"),
    title = "Convergence of finite differences",
)
```
``````

(demo-diffmats-cheb-julia)=
``````{dropdown} @demo-diffmats-cheb
Here is a $4\times 4$ Chebyshev differentiation matrix.

```{code-cell}
t, Dₓ = FNC.diffcheb(3, [-1, 1])
Dₓ
```

We again test the convergence rate.

```{code-cell}
f = x -> x + exp(sin(4 * x));
dfdx = x -> 1 + 4 * exp(sin(4 * x)) * cos(4 * x);
d2fdx2 = x -> 4 * exp(sin(4 * x)) * (4 * cos(4 * x)^2 - 4 * sin(4 * x));
```

```{code-cell}
n = 5:5:70
err1 = zeros(size(n))
err2 = zeros(size(n))
for (k, n) in enumerate(n)
    t, Dₓ, Dₓₓ = FNC.diffcheb(n, [-1, 1])
    y = f.(t)
    err1[k] = norm(dfdx.(t) - Dₓ * y, Inf)
    err2[k] = norm(d2fdx2.(t) - Dₓₓ * y, Inf)
end
```

Since we expect a spectral convergence rate, we use a semi-log plot for the error.

```{code-cell}
plot(
    n,
    [err1 err2],
    m = :o,
    label = [L"f'" L"f''"],
    xaxis = (L"n"),
    yaxis = (:log10, "max error"),
    title = "Convergence of Chebyshev derivatives",
)
```
``````

### 10.4 @section-bvp-linear
(demo-linear-solve-julia)=
``````{dropdown} @demo-linear-solve
We solve linear BVP  

$$ u'' - (\cos x) u' + (\sin x) u = 0, \quad u(0)=1, \; u\left(\frac{3\pi}{2}\right)=\frac{1}{e}. $$ 

Its exact solution is known:

```{code-cell}
exact = x -> exp(sin(x));
```

The problem is presented above in our standard form, so we can identify the coefficient functions in the ODE. Each should be coded as a function.

```{code-cell}
p = x -> -cos(x);
q = sin;
r = x -> 0;      # function, not value 
```

We solve the BVP and compare the result to the exact solution.

```{code-cell}
x, u = FNC.bvplin(p, q, r, [0, 3π / 2], 1, exp(-1), 30);
```

```{code-cell}
plot(exact, 0, 3π / 2, layout = (2, 1), label = "exact")
scatter!(
    x,
    u,
    m = :o,
    subplot = 1,
    label = "numerical",
    yaxis = ("solution"),
    title = "Solution of a linear BVP",
)

plot!(x, exact.(x) - u, subplot = 2, xaxis = L"x", yaxis = ("error"))
```
``````

(demo-linear-converge-julia)=
``````{dropdown} @demo-linear-converge
The BVP is
  
$$
u'' - \lambda^2 u = \lambda^2, \quad  u(0)=-1, \; u(1)=0,
$$

with exact solution $\sinh(\lambda x)/\sinh(\lambda) - 1$.

```{code-cell}
λ = 10
exact = x -> sinh(λ * x) / sinh(λ) - 1;
```

The following functions define the ODE.

```{code-cell}
p = x -> 0
q = x -> -λ^2
r = x -> λ^2;
```

We compare the computed solution to the exact one for increasing $n$.

```{code-cell}
n = 5 * [round(Int, 10^d) for d in 0:0.25:3]
err = zeros(size(n))
for (k, n) in enumerate(n)
    x, u = FNC.bvplin(p, q, r, [0, 1], -1, 0, n)
    err[k] = norm(exact.(x) - u, Inf)
end

data = (n = n[1:4:end], err = err[1:4:end])
pretty_table(data, header = ["n", "inf-norm error"])
```

Each factor of 10 in $n$ reduces error by a factor of 100, which is indicative of second-order convergence.

```{code-cell}
plot(
    n,
    err,
    m = :o,
    label = "observed",
    xaxis = (:log10, L"n"),
    yaxis = (:log10, "inf-norm error"),
    title = "Convergence for a linear BVP",
)
plot!(n, 0.25 * n .^ (-2), l = (:dash, :gray), label = "2nd order")
```
``````

### 10.5 @section-bvp-nonlinear
(demo-nonlinear-pendulum-julia)=
``````{dropdown} @demo-nonlinear-pendulum
Suppose a damped pendulum satisfies the nonlinear equation $\theta'' + 0.05\theta'+\sin \theta =0$. We want to start the pendulum at $\theta=2.5$ and give it the right initial velocity so that it reaches $\theta=-2$ at exactly $t=5$. This is a boundary-value problem with Dirichlet conditions $\theta(0)=2.5$ and $\theta(5)=-2$.

The first step is to define the function $\phi$ that equals $\theta''$.

```{code-cell}
ϕ = (t, θ, ω) -> -0.05 * ω - sin(θ);
```

Next, we define the boundary conditions.

```{code-cell}
g₁(u, du) = u - 2.5
g₂(u, du) = u + 2;
```

```{index} ! Julia; collect
```

::::{grid} 1 1 2 2

:::{grid-item}

The last ingredient is an initial estimate of the solution. Here we choose $n=100$ and a linear function between the endpoint values. 


:::

:::{card}

The `collect` function turns a range object into a true vector.

:::
::::

```{code-cell}
init = collect(range(2.5, -2, length = 101));
```

We find a solution with negative initial slope, i.e., the pendulum is initially pushed back toward equilibrium.

```{code-cell}
t, θ = FNC.bvp(ϕ, [0, 5], g₁, g₂, init)
plot(t, θ, xaxis = (L"t"), yaxis = (L"\theta(t)"), title = "Pendulum over [0,5]")
```

If we extend the time interval longer for the same boundary values, then the initial slope must adjust.

```{code-cell}
t, θ = FNC.bvp(ϕ, [0, 8], g₁, g₂, init)
plot(t, θ, xaxis = (L"t"), yaxis = (L"\theta(t)"), title = "Pendulum over [0,8]")
```

This time, the pendulum is initially pushed toward the unstable equilibrium in the upright vertical position before gravity pulls it back down.
``````

(demo-nonlinear-mems-julia)=
``````{dropdown} @demo-nonlinear-mems
We look for a solution to the parameterized membrane deflection problem from {numref}`Example {number} <example-tpbvp-mems>`,

$$
w''+ \frac{1}{r}w'= \frac{\lambda}{w^2},\quad w'(0)=0,\; w(1)=1.
$$ 

Here is the problem definition. We use a truncated domain to avoid division by zero at $r=0$.

```{code-cell}
domain = [eps(), 1]
λ = 0.5
ϕ = (r, w, dwdr) -> λ / w^2 - dwdr / r
g₁(w, dw) = dw
g₂(w, dw) = w - 1;
```

First we try a constant function as the initialization.

```{code-cell}
init = ones(301)
r, w₁ = FNC.bvp(ϕ, domain, g₁, g₂, init)

plot(r, w₁, xaxis = (L"r"), yaxis = (L"w(r)"), title = "Solution of the membrane problem")
```

It's not necessary that the initialization satisfy the boundary conditions. In fact, by choosing a different constant function as the initial guess, we arrive at another valid solution.

```{code-cell}
init = 0.5 * ones(301)
r, w₂ = FNC.bvp(ϕ, domain, g₁, g₂, init)
plot!(r, w₂, title = "Two solutions of the membrane problem")
```
``````

(demo-nonlinear-allencahn-julia)=
``````{dropdown} @demo-nonlinear-allencahn
We solve the stationary **Allen–Cahn equation**,
  
$$
\epsilon u'' = u^3-u, \quad 0 \le x \le 1, \quad u'(0)=0, \; u(1)=1.
$$

```{code-cell}
ϕ = (x, u, dudx) -> (u^3 - u) / ϵ;
g₁(u, du) = du
g₂(u, du) = u - 1;
```

Finding a solution is easy at larger values of $\epsilon$.

```{code-cell}
ϵ = 0.02
init = collect(range(-1, 1, length = 141))
x, u₁ = FNC.bvp(ϕ, [0, 1], g₁, g₂, init)

plot(
    x,
    u₁,
    label = L"\epsilon = 0.02",
    leg = :bottomright,
    xaxis = (L"x"),
    yaxis = (L"u(x)"),
    title = "Allen–Cahn solution",
)
```

However, finding a good initialization is not trivial for smaller values of $\epsilon$. Note below that the iteration stops without converging to a solution.

```{code-cell}
ϵ = 0.002;
x, z = FNC.bvp(ϕ, [0, 1], g₁, g₂, init);
```

The iteration succeeds if we use the first solution instead as the initialization here.

```{code-cell}
x, u₂ = FNC.bvp(ϕ, [0, 1], g₁, g₂, u₁)
plot!(x, u₂, label = L"\epsilon = 0.002")
```

In this case we can continue further.

```{code-cell}
ϵ = 0.0005
x, u₃ = FNC.bvp(ϕ, [0, 1], g₁, g₂, u₂)
plot!(x, u₃, label = L"\epsilon = 0.0005")
```
``````

### 10.6 @section-bvp-galerkin
(demo-galerkin-fem-julia)=
``````{dropdown} @demo-galerkin-fem
We solve the equation

$$
  -(x^2u')' + 4 y = \sin(\pi x), \qquad u(0)=u(1)=0,
$$

in which

$$
  c(x) = x^2, \qquad s(x) = 4, \qquad f(x)=\sin(\pi x).
$$

Here are the coefficient function definitions. Even though $s$ is a constant, it has to be defined as a function for {numref}`Function {number} <function-fem>` to use it.

```{code-cell}
c = x -> x^2;
q = x -> 4;
f = x -> sin(π * x);
```

```{code-cell}
x, u = FNC.fem(c, q, f, 0, 1, 50)
plot(
    x,
    u,
    label = "",
    xaxis = (L"x"),
    yaxis = (L"u"),
    title = "Solution by finite elements",
)
```
``````

