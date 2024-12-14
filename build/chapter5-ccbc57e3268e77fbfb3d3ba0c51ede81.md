---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---
```{code-cell}
:tags: [remove-cell]
import Pkg; Pkg.activate("/Users/driscoll/Documents/GitHub/fnc")
using FundamentalsNumericalComputation
FNC.init_format()
```

(demo-interpolation-global-julia)=
``````{dropdown} Trouble in polynomial interpolation
Here are some points that we could consider to be observations of an unknown function on $[-1,1]$.

```{code-cell}
n = 5
t = range(-1, 1, length=n + 1)
y = @. t^2 + t + 0.05 * sin(20 * t)

scatter(t, y, label="data", leg=:top)
```

```{index} ! Julia; fit
```

The polynomial interpolant, as computed using `fit`, looks very sensible. It's the kind of function you'd take home to meet your parents.

```{code-cell}
p = Polynomials.fit(t, y, n)     # interpolating polynomial
plot!(p, -1, 1, label="interpolant")
```

But now consider a different set of points generated in almost exactly the same way.

```{code-cell}
n = 18
t = range(-1, 1, length=n + 1)
y = @. t^2 + t + 0.05 * sin(20 * t)

scatter(t, y, label="data", leg=:top)
```

The points themselves are unremarkable. But take a look at what happens to the polynomial interpolant.

```{code-cell}
p = Polynomials.fit(t, y, n)
x = range(-1, 1, length=1000)    # use a lot of points
plot!(x, p.(x), label="interpolant")
```

Surely there must be functions that are more intuitively representative of those points!
``````

(demo-interpolation-pwise-julia)=
``````{dropdown} Piecewise polynomial interpolation
Let us recall the data from {numref}`Demo %s <demo-interpolation-global>`.

```{code-cell}
n = 12
t = range(-1, 1, length=n + 1)
y = @. t^2 + t + 0.5 * sin(20 * t)

scatter(t, y, label="data", leg=:top)
```

Here is an interpolant that is linear between each consecutive pair of nodes, using `plinterp` from {numref}`section-localapprox-pwlin`.

```{code-cell}
p = FNC.plinterp(t, y)
plot!(p, -1, 1, label="piecewise linear")
```

```{index} ! Julia; Spline1D
```

We may prefer a smoother interpolant that is piecewise cubic, generated using `Spline1D` from the `Dierckx` package.

```{code-cell}
p = Spline1D(t, y)
plot!(x -> p(x), -1, 1, label="piecewise cubic")
```
``````

(demo-interp-cond-julia)=
``````{dropdown} Conditioning of interpolation
In {numref}`Demo %s <demo-interpolation-global>` and {numref}`Demo %s <demo-interpolation-pwise>` we saw a big difference between polynomial interpolation and piecewise polynomial interpolation of some arbitrarily chosen data. The same effects can be seen clearly in the cardinal functions, which are closely tied to the condition numbers.

```{code-cell}
n = 18
t = range(-1, stop=1, length=n + 1)
y = [zeros(9); 1; zeros(n - 9)];  # data for 10th cardinal function

scatter(t, y, label="data")
```

```{code-cell}
ϕ = Spline1D(t, y)
plot!(x -> ϕ(x), -1, 1, label="spline",
    xlabel=L"x", ylabel=L"\phi(x)",
    title="Piecewise cubic cardinal function")
```

The piecewise cubic cardinal function is nowhere greater than one in absolute value. This happens to be true for all the cardinal functions, ensuring a good condition number for any interpolation with these functions. But the story for global polynomials is very different.

```{code-cell}
scatter(t, y, label="data")

ϕ = Polynomials.fit(t, y, n)
plot!(x -> ϕ(x), -1, 1, label="polynomial",
    xlabel=L"x", ylabel=L"\phi(x)", legend=:top,
    title="Polynomial cardinal function")
```

From the figure we can see that the condition number for polynomial interpolation on these nodes is at least 500.
``````

(function-hatfun-julia)=
``````{dropdown} Hat function
```{literalinclude} ../julia/package/src/chapter05.jl
:filename: hatfun.jl
:start-line: 0
:end-line: 19
:language: julia
:linenos: true
```
``````

(demo-pwlin-hat-julia)=
``````{dropdown} A look at hat functions
Let's define a set of four nodes (i.e., $n=3$ in our formulas).

```{index} ! Julia; annotate!
```

```{code-cell}
t = [0, 0.55, 0.7, 1]
```

::::{grid} 1 1 2 2
We plot the hat functions $H_0,\ldots,H_3$.
:::{card}
Use `annotate!` to add text to a plot.
:::
::::

```{code-cell}
plt = plot(layout=(4, 1), legend=:top,
    xlabel=L"x", ylims=[-0.1, 1.1], ytick=[])
for k in 0:3
    Hₖ = FNC.hatfun(t, k)
    plot!(Hₖ, 0, 1, subplot=k + 1)
    scatter!(t, Hₖ.(t), m=3, subplot=k + 1)
    annotate!(t[k+1], 0.25, text(latexstring("H_$k"), 10), subplot=k + 1)
end
plt
```
``````

(demo-pwlin-usage-julia)=
``````{dropdown} Using piecewise linear interpolation
We generate a piecewise linear interpolant of $f(x)=e^{\sin 7x}$.

```{code-cell}
f = x -> exp(sin(7 * x))

plot(f, 0, 1, label="function", xlabel=L"x", ylabel=L"y")
```

First we sample the function to create the data.

```{code-cell}
t = [0, 0.075, 0.25, 0.55, 0.7, 1]    # nodes
y = f.(t)                             # function values

scatter!(t, y, label="values at nodes")
```

Now we create a callable function that will evaluate the piecewise linear interpolant at any $x$, and then plot it.

```{code-cell}
p = FNC.plinterp(t, y)
plot!(p, 0, 1, label="interpolant", title="PL interpolation")
```
``````

(function-plinterp-julia)=
``````{dropdown} Piecewise linear interpolation
```{literalinclude} ../julia/package/src/chapter05.jl
:filename: plinterp.jl
:start-line: 28
:end-line: 37
:language: julia
:linenos: true
```
``````

