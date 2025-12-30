---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.12
numbering:
  headings: false
---
# Chapter 4

## Functions

(function-newton-julia)=
``````{dropdown} Newton's method
:open:
```{literalinclude} FNCFunctions/src/chapter04.jl
:filename: newton.jl
:start-after: # begin newton
:end-before: # end newton
:language: julia
:linenos: true
```
```{admonition} About the code
:class: dropdown
{numref}`Function {number} <function-newton>` accepts *keyword arguments*. In the function declaration, these follow the semicolon, and when the function is called, they may be supplied as `keyword=value` in the argument list. Here, these arguments are also given default values by the assignments within the declaration. This arrangement is useful when there are multiple optional arguments, because the ordering of them doesn't matter.

The `break` statement, seen here in line 25, causes an immediate exit from the innermost loop in which it is called. It is often used as a safety valve to escape an iteration that may not be able to terminate otherwise.
```
``````

(function-secant-julia)=
``````{dropdown} Secant method
:open:
```{literalinclude} FNCFunctions/src/chapter04.jl
:filename: secant.jl
:start-after: # begin secant
:end-before: # end secant
:language: julia
:linenos: true
```
```{admonition} About the code
:class: dropdown
Because we want to observe the convergence of the method, {numref}`Function {number} <function-secant>` stores and returns the entire sequence of root estimates. However, only the most recent two are needed by the iterative formula. This is demonstrated by the use of `y₁` and `y₂` for the two most recent values of $f$.
```
``````

(function-newtonsys-julia)=
``````{dropdown} Newton's method for systems
:open:
```{literalinclude} FNCFunctions/src/chapter04.jl
:filename: newtonsys.jl
:start-after: # begin newtonsys
:end-before: # end newtonsys
:language: julia
:linenos: true
```
````{admonition} About the code
:class: dropdown
The output of {numref}`Function {number} <function-newtonsys>` is a vector of vectors representing the entire history of root estimates. Since these should be in floating point, the starting value is converted with `float` before the iteration starts.
````
``````

(function-fdjac-julia)=
``````{dropdown} Finite differences for Jacobian
:open:
```{literalinclude} FNCFunctions/src/chapter04.jl
:filename: fdjac.jl
:start-after: # begin fdjac
:end-before: # end fdjac
:language: julia
:linenos: true
```
::::{admonition} About the code
:class: dropdown
{numref}`Function {number} <function-fdjac>` is written to accept the case where $\mathbf{f}$ maps $n$ variables to $m$ values with $m\neq n$, in anticipation of {numref}`section-nonlineqn-nlsq`.

Note that a default value is given for the third argument `y₀`, and it refers to earlier arguments in the list. The reason is that in some contexts, the caller of `fdjac` may have already computed `y₀` and can supply it without computational cost, while in other contexts, it must be computed fresh. The configuration here adapts to either situation.
::::
``````

(function-levenberg-julia)=
``````{dropdown} Levenberg's method
:open:
```{literalinclude} FNCFunctions/src/chapter04.jl
:filename: levenberg.jl
:start-after: # begin levenberg
:end-before: # end levenberg
:language: julia
:linenos: true
```
``````

## Examples

```{code-cell}
:tags: remove-cell
include("FNC_init.jl")
```
### 4.1 @section-nonlineqn-rootproblem
(demo-rootproblem-bessel-julia)=
``````{dropdown} @demo-rootproblem-bessel
:open:

```{code-cell}
using SpecialFunctions
J₃(x) = besselj(3, x)
fig, ax, plt = lines(0..20, J₃;
    axis=(title="Bessel function", xlabel=L"x", ylabel=L"J_3(x)")
)
```
From the graph we see roots near 6, 10, 13, 16, and 19. We use `nlsolve` from the `NLsolve` package to find these roots accurately. It uses vector variables, so we have to code accordingly.
```{tip}
:class: dropdown
Type `\omega` followed by <kbd>Tab</kbd> to get the character `ω`.
The argument `ftol=1e-14` below is called a **keyword argument**. Here it sets a goal for the maximum value of $|f(x)|$.
```

```{code-cell}
using NLsolve
ω = []
for guess = [6., 10. ,13., 16., 19.]
    s = nlsolve(x -> J₃(x[1]), [guess], ftol=1e-14)
    append!(ω, s.zero)
end
```

```{code-cell}
y = J₃.(ω)
pretty_table( [ω y];
    column_labels=["root estimate", "function value"], backend=:html) 
```

```{code-cell}
scatter!(ω, y, color=:red)
fig
```

If instead we seek values at which $J_3(x)=0.2$, then we must find roots of the function $J_3(x)-0.2$.

```{code-cell}
r = []
for guess = [3., 6., 10., 13.]
    f(x) = J₃(x[1]) - 0.2
    s = nlsolve(f, [guess], ftol=1e-14)
    append!(r, s.zero)
end
scatter!(r, J₃.(r), color=:black)
fig
```
``````

(demo-roots-cond-julia)=
``````{dropdown} @demo-roots-cond
:open:
Consider first the function

```{code-cell}
f(x) = (x - 1) * (x - 2);
```

:::{index} ! Julia; splatting
:::


At the root $r=1$, we have $f'(r)=-1$. If the values of $f$ were perturbed at every point by a small amount of noise, we can imagine finding the root of the function drawn with a thick ribbon, giving a range of potential roots.
```{tip}
:class: dropdown
The syntax `interval...` is called **splatting** and means to insert all the individual elements of `interval` as a sequence.
```


```{code-cell}
x = range(0.8, 1.2, 400)
y = f.(x)
fig = Figure()
ax = Axis(fig[1, 1];
    xlabel=L"x", ylabel=L"f(x)", title="Well-conditioned root", aspect=DataAspect())
band!(x, y .- 0.03, y .+ 0.03, color=:lightblue, alpha=0.5)
lines!(x, y) 
scatter!([1], [0])
ylims!(ax, -0.2, 0.2)
fig
```

The possible values for a perturbed root all lie within the interval where the ribbon intersects the $x$-axis. The width of that zone is about the same as the vertical thickness of the ribbon.

By contrast, consider the function

```{code-cell}
f(x) = (x - 1) * (x - 1.01);
```

Now $f'(1)=-0.01$, and the graph of $f$ will be much shallower near $x=1$. Look at the effect this has on our thick rendering:

```{code-cell}
x = range(0.8, 1.2, 400)
y = f.(x)
fig = Figure()
ax = Axis(fig[1, 1];
    xlabel=L"x", ylabel=L"f(x)", title="Ill-conditioned root", aspect=DataAspect())
band!(x, y .- 0.03, y .+ 0.03, color=:lightblue, alpha=0.5)
lines!(x, y) 
scatter!([1], [0])
ylims!(ax, -0.2, 0.2)
fig
```

The vertical displacements in this picture are exactly the same as before. But the potential _horizontal_ displacement of the root is much wider. In fact, if we perturb the function entirely upward by the amount drawn here, the root disappears!

``````

### 4.2 @section-nonlineqn-fixed-point

(demo-fp-spiral-julia)=
``````{dropdown} @demo-fp-spiral
:open:
Let's convert the roots of a quadratic polynomial $f(x)$ to a fixed point problem.

```{code-cell}
using Polynomials
p = Polynomial([3.5, -4,1])
r = roots(p)
rmin, rmax = extrema(r)
@show rmin, rmax;
```

We define $g(x)=x-p(x)$. 

```{code-cell}
g(x) = x - p(x)
```

Intersections of $y=g(x)$ with the line $y=x$ are fixed points of $g$ and thus roots of $f$. (Only one is shown in the chosen plot range.)

```{code-cell}
fig = Figure()
ax = Axis(fig[1, 1];
    xlabel=L"x", ylabel=L"y", title="Finding a fixed point", 
    autolimitaspect=1, xgridvisible=false, ygridvisible=false)
lines!(2..2.8, g; label=L"y=g(x)")
lines!(2..2.8, identity; label=L"y=x", color=:black, linestyle=:dash)
ylims!(ax, 2.4, 2.8)
axislegend(ax, position=:rb)
fig
```

If we evaluate $g(2.1)$, we get a value of almost 2.6, so this is not a fixed point.

```{code-cell}
x = 2.1;
y = g(x)
```

However, $y=g(x)$ is considerably closer to the fixed point at around 2.7 than $x$ is. Suppose then that we adopt $y$ as our new $x$ value. Changing the $x$ coordinate in this way is the same as following a horizontal line over to the graph of $y=x$.

```{code-cell}
arrow(A, B, color) = arrows2d!(A, B; color, argmode=:endpoints, shaftwidth=2)
arrow(Point(x, y), Point(y, y), :purple)
fig
```

Now we can compute a new value for $y$. We leave $x$ alone here, so we travel along a vertical line to the graph of $g$.

```{code-cell}
x = y;  y = g(x)
arrow(Point(x, x), Point(x, y), :orange)
fig
```

You see that we are in a position to repeat these steps as often as we like. Let's apply them a few times and see the result.

```{code-cell}
for k = 1:5
    arrow(Point(x, y), Point(y, y), :purple)
    x = y       # y becomes the new x
    y = g(x)    # g(x) becomes the new y
    arrow(Point(x, x), Point(x, y), :orange)
end
fig
```

The process spirals in beautifully toward the fixed point we seek. Our last estimate has almost 4 accurate digits.

```{code-cell} 
abs(y - rmax) / rmax
```

Now let's try to find the other fixed point $\approx 1.29$ in the same way. We'll use 1.3 as a starting approximation.

```{code-cell}
:tags: [hide-input]
fig = Figure()
ax = Axis(fig[1, 1];
    xlabel=L"x", ylabel=L"y", title="Finding a fixed point", aspect=DataAspect())
lines!(1..2, g; label=L"y=g(x)")
lines!(1..2, identity; label=L"y=x", color=:black, linestyle=:dash)
axislegend(ax, position=:rb)

x = 1.32; y = g(x);
for k = 1:4
    arrow(Point(x, y), Point(y, y), :purple)
    x = y       # y --> new x
    y = g(x)    # g(x) --> new y
    arrow(Point(x, x), Point(x, y), :orange)
end
fig
```

This time, the iteration is pushing us _away from_ the correct answer.
``````

(demo-fp-converge-julia)=
``````{dropdown} @demo-fp-converge
:open:
We revisit @demo-fp-spiral and investigate the observed convergence more closely. Recall that above we calculated $g'(p)\approx-0.42$ at the convergent fixed point.

```{code-cell}
p = Polynomial([3.5, -4, 1])
r = roots(p)
rmin, rmax = extrema(r)
@show rmin, rmax;
```

Here is the fixed-point iteration. This time we keep track of the whole sequence of approximations.

:::{index} Julia; push!
:::

```{code-cell}
g(x) = x - p(x)
x = [2.1]
for k = 1:12
    push!(x, g(x[k]))
end
x
```

It's illuminating to construct and plot the sequence of errors.

```{code-cell}
err = @. abs(x - rmax)
scatterlines(0:12, err;
    axis=(xlabel="iteration number",  ylabel="error", yscale=log10,
    title="Convergence of fixed-point iteration"))
```

It's quite clear that the convergence quickly settles into a linear rate. We could estimate this rate by doing a least-squares fit to a straight line. Keep in mind that the values for small $k$ should be left out of the computation, as they don't represent the linear trend.

```{code-cell}
y = log.(err[5:12])
p = Polynomials.fit(5:12, y, 1)
```

We can exponentiate the slope to get the convergence constant $\sigma$.

```{code-cell}
σ = exp(p.coeffs[2])
```

The error should therefore decrease by a factor of $\sigma$ at each iteration. We can check this easily from the observed data.

```{code-cell}
[err[i+1] / err[i] for i in 8:11]
```

The methods for finding $\sigma$ agree well.
``````

### 4.3 @section-nonlineqn-newton

(demo-newton-line-julia)=
``````{dropdown} @demo-newton-line
:open:

Suppose we want to find a root of the function $f(x)=x e^x - 2$.

```{code-cell}
f(x) = x * exp(x) - 2
fig, ax, plt = lines(0..1.5, f; label="function",
    axis=(xgridvisible=false, xlabel=L"x",  ylabel=L"y"))
ylims!(ax, -2, 4)
fig
```

From the graph, it is clear that there is a root near $x=1$. So we call that our initial guess, $x_1$.

```{code-cell}
x₁ = 1
y₁ = f(x₁)
scatter!([x₁], [y₁], label="initial point")
text!(x₁ + 0.02, y₁, text=L"(x_1, f(x_1))", align = (:left, :center))
fig
```

Next, we can compute the tangent line at the point $\bigl(x_1,f(x_1)\bigr)$, using the derivative.

```{code-cell}
df_dx(x) = exp(x) * (x + 1)
m₁ = df_dx(x₁)
tangent = x -> y₁ + m₁ * (x - x₁)

lines!(0..1.5, tangent; linestyle=:dash, label="tangent line")
axislegend(ax, position=:lt)
ax.title = "Tangent-line approximation"
fig
```

In lieu of finding the root of $f$ itself, we settle for finding the root of the tangent line approximation, which is trivial. Call this $x_2$, our next approximation to the root.

```{code-cell}
@show x₂ = x₁ - y₁ / m₁
scatter!([x₂], [0])
text!(x₂ + 0.02, 0, text=L"(x_2, 0)", align = (:left, :center))
ax.title = "First iteration"
fig
```

```{code-cell}
y₂ = f(x₂)
```

The residual (i.e., value of $f$) is smaller than before, but not zero. So we repeat the process with a new tangent line based on the latest point on the curve.

```{code-cell}
:tags: [hide-input]
fig, ax, plt = lines(0.82..0.87, f; label="function",
    axis=(xgridvisible=false, xlabel=L"x",  ylabel=L"y", title="Second iteration"))

scatter!([x₂], [y₂], label="starting point")
text!(x₂ + 0.002, y₂, text=L"(x_2, f(x_2))", align = (:left, :center))

m₂ = df_dx(x₂)
tangent = x -> y₂ + m₂ * (x - x₂)
lines!(0.82..0.87, tangent; linestyle=:dash, label="tangent line")

@show x₃ = x₂ - y₂ / m₂
scatter!([x₃], [0], label="tangent root")
text!(x₃ - 0.002, 0, text=L"(x_3, 0)", align = (:right, :center))
axislegend(ax, position=:lt)
fig
```

```{code-cell}
y₃ = f(x₃)
```

Judging by the residual, we appear to be getting closer to the true root each time.
``````

(demo-newton-converge-julia)=
``````{dropdown} @demo-newton-converge
:open:
We again look at finding a solution of $x e^x=2$ near $x=1$. To apply Newton's method, we need to calculate values of both the residual function $f$ and its derivative.

```{code-cell}
f(x) = x * exp(x) - 2;
df_dx(x) = exp(x) * (x + 1);
```

We don't know the exact root, so we use `nlsolve` to determine a proxy for it.

```{code-cell}
using NLsolve
r = nlsolve(x -> f(x[1]), [1.0]).zero
```

We use $x_1=1$ as a starting guess and apply the iteration in a loop, storing the sequence of iterates in a vector.

```{code-cell}
x = [1; zeros(4)]
for k = 1:4
    x[k+1] = x[k] - f(x[k]) / df_dx(x[k])
end
x
```

Here is the sequence of errors.

```{code-cell}
ϵ = @. x - r
```

Because the error reaches machine epsilon so rapidly, we're going to use extended precision to allow us to take a few more iterations. We'll take the last iteration as the most accurate root estimate.
```{tip}
:class: dropdown
A `BigFloat` uses 256 bits of precision, rather than 53 in `Float64`. But arithmetic is done by software emulation and is much slower.
```

```{code-cell}
x = [BigFloat(1); zeros(7)]
for k = 1:7
    x[k+1] = x[k] - f(x[k]) / df_dx(x[k])
end
r = x[end]
```

```{code-cell}
ϵ = @. Float64(x[1:end-1] - r)
```

The exponents in the scientific notation definitely suggest a squaring sequence. We can check the evolution of the ratio in {eq}`quadratictest`.

```{code-cell}
logerr = @. log10(abs(ϵ))
ratios = [NaN; [logerr[i+1] / logerr[i] for i in 1:length(logerr)-1]]
pretty_table( (iter=1:7, ϵ, logerr, ratios);
    column_labels=["iteration", "error", "log error", "ratio"], backend=:html )
```

The clear convergence to 2 above constitutes good evidence of quadratic convergence.
``````

(demo-newton-usage-julia)=
``````{dropdown} @demo-newton-usage
:open:
```{index} ! Julia; enumerate
```

Suppose we want to evaluate the inverse of the function $h(x)=e^x-x$. This means solving $y=e^x-x$ for $x$ when $y$ is given, which has no elementary form. If a value of $y$ is given numerically, though, we simply have a rootfinding problem for $f(x)=e^x-x-y$.
```{tip}
:class: dropdown
The `enumerate` function produces a pair of values for each iteration: a positional index and the corresponding contents.
```

```{code-cell}
g(x) = exp(x) - x
dg_dx(x) = exp(x) - 1
y = range(g(0), g(2), 200)
x = zeros(length(y))
for (i, y) in enumerate(y)
    f(x) = g(x) - y
    df_dx(x) = dg_dx(x)
    r = FNC.newton(f, df_dx, y)
    x[i] = r[end]
end

fig, ax, plt = lines(0..2, g; label=L"g(x)",
    axis=(xlabel=L"x",  ylabel=L"y", title="Function and its inverse", aspect=DataAspect()))
lines!(y, x, label=L"g^{-1}(y)")
lines!(0..maximum(y), identity, linestyle=:dash, color=:black)
axislegend(ax, position=:rb)
fig
```
``````

### 4.4 @section-nonlineqn-secant
(demo-secant-line-julia)=
``````{dropdown} @demo-secant-line
:open:


We return to finding a root of the equation $x e^x=2$.

```{code-cell}
f(x) = x * exp(x) - 2;
fig, ax, plt = lines(0.25..1.25, f; label="function",
    axis=(xlabel=L"x",  ylabel=L"y", title="Two initial values"))
```

From the graph, it's clear that there is a root near $x=1$. To be more precise, there is a root in the interval $[0.5,1]$. So let us take the endpoints of that interval as _two_ initial approximations.

```{code-cell}
x₁ = 1;
y₁ = f(x₁);
x₂ = 0.5;
y₂ = f(x₂);
scatter!([x₁, x₂], [y₁, y₂]; label="initial points")
fig
```

Instead of constructing the tangent line by evaluating the derivative, we can construct a linear model function by drawing the line between the two points $\bigl(x_1,f(x_1)\bigr)$ and $\bigl(x_2,f(x_2)\bigr)$. This is called a _secant line_.

```{code-cell}
m₂ = (y₂ - y₁) / (x₂ - x₁)
secant = x -> y₂ + m₂ * (x - x₂)
lines!(0.25..1.25, secant; label="secant line", linestyle=:dash, color=:black)
axislegend(ax, position=:lt)
fig
```

As before, the next root estimate in the iteration is the root of this linear model.

```{code-cell}
x₃ = x₂ - y₂ / m₂
@show y₃ = f(x₃)
scatter!([x₃], [0])
text!(x₃ + 0.02, 0, text=L"(x_3, 0)", align = (:left, :center))
ax.title= "First iteration"
fig
```

For the next linear model, we use the line through the two most recent points. The next iterate is the root of that secant line, and so on.

```{code-cell}
m₃ = (y₃ - y₂) / (x₃ - x₂)
x₄ = x₃ - y₃ / m₃
```
``````

(demo-secant-converge-julia)=
``````{dropdown} @demo-secant-converge
:open:
We check the convergence of the secant method from @demo-secant-line. Again we will use extended precision to get a longer sequence than double precision allows.

```{code-cell}
f(x) = x * exp(x) - 2
x = FNC.secant(f, BigFloat(1), BigFloat(0.5), xtol=1e-80, ftol=1e-80);
```

We don't know the exact root, so we use the last value as a proxy.

```{code-cell}
r = x[end]
```

Here is the sequence of errors.

```{code-cell}
ϵ = @. Float64(r - x[1:end-2])
```

It's not easy to see the convergence rate by staring at these numbers. We can use {eq}`superlinear-rate` to try to expose the superlinear convergence rate.

```{code-cell}
logerr = @. log10(abs(ϵ))
ratios = [NaN; [logerr[i+1] / logerr[i] for i in 1:length(logerr)-1]]
pretty_table( (iter=eachindex(ϵ), ϵ, logerr, ratios);
    column_labels=["iteration", "error", "log error", "ratio"], backend=:html )
```

As expected, this settles in at around 1.618.
``````

(demo-secant-iqi-julia)=
``````{dropdown} @demo-secant-iqi
:open:
Here we look for a root of $x+\cos(10x)$ that is close to 1.

```{code-cell}
f(x) = x + cos(10 * x)
pts = range(0.5, 1.5, 500)

fig = Figure()
ax = Axis(fig[1, 1]; xgridvisible=false, xlabel=L"x", ylabel=L"y")
lines!(pts, f.(pts); label="function")
ylims!(ax, -0.1, 3)
fig
```
We choose three values to get the iteration started.

```{code-cell}
x = [0.8, 1.2, 1]
y = @. f(x)
scatter!(x, y, color=:black, label="initial points")
fig
```

If we were using forward interpolation, we would ask for the polynomial interpolant of $y$ as a function of $x$. But that parabola has no real roots.

```{code-cell}
using Polynomials
q = Polynomials.fit(x, y, 2)      # interpolating polynomial
lines!(pts, q.(pts); linestyle=:dash, label="interpolant")
fig
```

To do inverse interpolation, we swap the roles of $x$ and $y$ in the interpolation.

```{code-cell}
fig = Figure()
ax = Axis(fig[1, 1]; xgridvisible=false, xlabel=L"x", ylabel=L"y")
lines!(pts, f.(pts); label="function")
scatter!(x, y, color=:black, label="initial points")

vals = range(-0.3, 2.2, 500)
q = Polynomials.fit(y, x, 2)       # interpolating polynomial
lines!(q.(vals), vals; linestyle=:dash, label="inverse interpolant")
axislegend(ax, position=:rc)
fig
```

We seek the value of $x$ that makes $y$ zero. This means evaluating $q$ at zero.

```{code-cell}
q(0)
```

Let's restart the process with `BigFloat` numbers to get more precision in order to see the convergence.

```{code-cell}
x = BigFloat.([8, 12, 10]) / 10
y = @. f(x)

for k = 3:12
    q = Polynomials.fit(y[k-2:k], x[k-2:k], 2)
    push!(x, q(0))
    push!(y, f(x[k+1]))
end

println("residual = $(f(x[end]))")
```

As far as our current precision is concerned, we have an exact root.

```{code-cell}
r = x[end]
ϵ = @. Float64(abs(r - x[1:end-1]))
logerr = @. log10(abs(ϵ))
ratios = [NaN; [logerr[i+1] / logerr[i] for i in 1:length(logerr)-1]]
pretty_table( (iter=eachindex(ϵ), ϵ, logerr, ratios);
    column_labels=["iteration", "error", "log error", "ratio"], backend=:html )
```

The convergence is apparently superlinear at a rate of $\alpha=1.8$ or so.
``````

### 4.5 @section-nonlineqn-newtonsys

(demo-newtonsys-converge-julia)=
``````{dropdown} @demo-newtonsys-converge
:open:
A system of nonlinear equations is defined by its residual and Jacobian.
```{tip}
:class: dropdown
Be careful when coding a Jacobian all in one statement. Spaces separate columns, so `x[3]-1` is not the same as `x[3] - 1`.
```

```{code-cell}
function func(x)
    return [exp(x[2] - x[1]) - 2,
        x[1] * x[2] + x[3],
        x[2] * x[3] + x[1]^2 - x[2]
    ]
end;

function jac(x)
    return [
        -exp(x[2]-x[1]) exp(x[2]-x[1]) 0
        x[2] x[1] 1
        2x[1] x[3]-1 x[2]
    ]
end;
```

We will use a `BigFloat` starting value, and commensurately small stopping tolerances, in order to get a sequence long enough to measure convergence.

```{code-cell}
x₁ = BigFloat.([0, 0, 0])
ϵ = eps(BigFloat)
x = FNC.newtonsys(func, jac, x₁, xtol=ϵ, ftol=ϵ);
```

Let's compute the residual of the last result in order to check the quality.

```{code-cell}
r = x[end]
@show residual = norm(func(r));
```

We take the sequence of norms of errors, applying the log so that we can look at the exponents.

```{code-cell}
logerr = [Float64(log(norm(r - x[k]))) for k in 1:length(x)-1]
ratios = [NaN; [logerr[i+1] / logerr[i] for i in 1:length(logerr)-1]]
pretty_table( (iter=eachindex(logerr), logerr, ratios);
    column_labels=["iteration", "log error", "ratio"], backend=:html )
```

The ratio is neatly converging toward 2, which is expected for quadratic convergence.
``````

### 4.6 @section-nonlineqn-quasinewton

(demo-quasi-levenberg-julia)=
``````{dropdown} @demo-quasi-levenberg
:open:
To solve a nonlinear system, we need to code only the function defining the system, and not its Jacobian.

```{code-cell}
f(x) = 
    return [
        exp(x[2] - x[1]) - 2,
        x[1] * x[2] + x[3],
        x[2] * x[3] + x[1]^2 - x[2]
    ]
```

In all other respects usage is the same as for the `newtonsys` function.

```{code-cell}
x₁ = [0.0, 0.0, 0.0]
x = FNC.levenberg(f, x₁)
```

It's always a good idea to check the accuracy of the root, by measuring the residual (backward error).

```{code-cell}
r = x[end]
println("backward error = $(norm(f(r)))")
```

Looking at the convergence in norm, we find a convergence rate between linear and quadratic, like with the secant method.

```{code-cell}
logerr = [log(norm(r - x[k])) for k in 1:length(x)-1]
ratios = [NaN; [logerr[i+1] / logerr[i] for i in 1:length(logerr)-1]]
pretty_table((iter=eachindex(logerr), logerr, ratios);
    column_labels=["iteration", "log error", "ratio"], backend=:html )
```
``````

### 4.7 @section-nonlineqn-nlsq

(demo-nlsq-converge-julia)=
``````{dropdown} @demo-nlsq-converge
:open:
We will observe the convergence of {numref}`Function {number} <function-levenberg>` for different levels of the minimum least-squares residual. We start with a function mapping from $\real^2$ into $\real^3$, and a point that will be near the optimum.

```{code-cell}
g(x) = [sin(x[1] + x[2]), cos(x[1] - x[2]), exp(x[1] - x[2])]
p = [1, 1];
```

```{index} ! Julia; @sprintf
```

The function $\mathbf{g}(\mathbf{x}) - \mathbf{g}(\mathbf{p})$ obviously has a zero residual at $\mathbf{p}$. We'll make different perturbations of that function in order to create nonzero residuals.
```{tip}
:class: dropdown
`@sprintf` is a way to format numerical values as strings, patterned after the C function `printf`.
```

```{code-cell}
using Printf
fig = Figure()
ax = Axis(fig[1, 1]; title="Convergence of Gauss–Newton",
    xlabel="iteration", ylabel="error", yscale=log10)
for R in [1e-3, 1e-2, 1e-1]
    # Define the perturbed function.
    f(x) = g(x) - g(p) + R * normalize([-1, 1, -1])
    x = FNC.levenberg(f, [0, 0])
    r = x[end]
    err = [norm(x - r) for x in x[1:end-1]]
    normres = norm(f(r))
    lines!(err, label=@sprintf("R=%.2g", normres))
end
axislegend(ax, position=:rt)
fig
```

In the least perturbed case, where the minimized residual is less than $10^{-3}$, the convergence is plausibly quadratic. At the next level up, the convergence starts similarly but suddenly stagnates for a long time. In the most perturbed case, the quadratic phase is nearly gone and the overall shape looks linear.
``````

(demo-nlsq-MM-julia)=
``````{dropdown} @demo-nlsq-MM
:open:
```{code-cell}
m = 25;
s = range(0.05, 6, length=m)
ŵ = @. 2 * s / (0.5 + s)                      # exactly on the curve
w = @. ŵ + 0.15 * cos(2 * exp(s / 16) * s);     # smooth noise added
```

```{code-cell}
fig, ax, plt = scatter(s, w; color=:black, label="noisy data",
    axis=(xlabel="s", ylabel="v"))
lines!(s, ŵ; linestyle=:dash, color=:black, label="perfect data")
leg = axislegend(ax, position=:rb)
fig
```

```{index} ! Julia; destructuring
```

The idea is to pretend that we know nothing of the origins of this data and use nonlinear least squares to find the parameters in the theoretical model function $v(s)$. In {eq}`nlsq-misfit`, the $s$ variable plays the role of $t$, and $v$ plays the role of $g$.
```{tip}
:class: dropdown
Putting comma-separated values on the left of an assignment will **destructure** the right-hand side, drawing individual assignments from entries of a vector, for example.
```

```{code-cell}
function misfit(x)
    V, Km = x   # rename components for clarity
    return @. V * s / (Km + s) - w
end
```

In the Jacobian the derivatives are with respect to the parameters in $\mathbf{x}$.

```{code-cell}
function misfitjac(x)
    V, Km = x   # rename components for clarity
    J = zeros(m, 2)
    J[:, 1] = @. s / (Km + s)              # dw/dV
    J[:, 2] = @. -V * s / (Km + s)^2         # dw/d_Km
    return J
end
```

```{code-cell}
x₁ = [1, 0.75]
x = FNC.newtonsys(misfit, misfitjac, x₁)

@show V, Km = x[end]  # final values
```

The final values are reasonably close to the values $V=2$, $K_m=0.5$ that we used to generate the noise-free data. Graphically, the model looks close to the original data.

```{code-cell}
model(s) = V * s / (Km + s)
lines!(0..6, model; label="nonlinear fit")
delete!(leg)
leg = axislegend(ax, position=:rb)
fig
```

For this particular model, we also have the option of linearizing the fit process. Rewrite the model as 

```{math}
:enumerated: false
\frac{1}{w} = \frac{\alpha}{s} + \beta = \alpha \cdot s^{-1} + \beta
```

for the new fitting parameters $\alpha=K_m/V$ and $\beta=1/V$. This corresponds to the misfit function whose entries are

$$f_i([\alpha,\beta]) = \left(\alpha \cdot \frac{1}{s_i} + \beta\right) - \frac{1}{w_i}$$
for $i=1,\ldots,m$. Although this misfit is nonlinear in $s$ and $w$, it's linear in the unknown parameters $\alpha$ and $\beta$. This lets us pose and solve it as a linear least-squares problem.

```{code-cell}
A = [s .^ (-1) s .^ 0]
u = 1 ./ w
α, β = A \ u
```

The two fits are different because they do not optimize the same quantities.

```{code-cell}
linmodel(x) = 1 / (β + α / x)
lines!(0..6, linmodel; label="linearized fit")
delete!(leg)
axislegend(ax, position=:rb)
fig
```

The truly nonlinear fit is clearly better in this case. It optimizes a residual for the original measured quantity rather than a transformed one we picked for algorithmic convenience.
``````
