---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{code-cell}
include("../../../julia/FNC_init.jl")
```
[**Demo %s**](#demo-secant-iqi)

Here we look for a root of $x+\cos(10x)$ that is close to 1.

```{code-cell}
f(x) = x + cos(10 * x)
interval = [0.5, 1.5]

plot(f, interval..., label="Function", legend=:bottomright,
    grid=:y, ylim=[-0.1, 3], xlabel=L"x", ylabel=L"y")
```

We choose three values to get the iteration started.

```{code-cell}
x = [0.8, 1.2, 1]
y = @. f(x)
scatter!(x, y, label="initial points")
```

If we were using forward interpolation, we would ask for the polynomial interpolant of $y$ as a function of $x$. But that parabola has no real roots.

```{code-cell}
using Polynomials
q = Polynomials.fit(x, y, 2)      # interpolating polynomial
plot!(x -> q(x), interval..., l=:dash, label="interpolant")
```

To do inverse interpolation, we swap the roles of $x$ and $y$ in the interpolation.
:::
```{tip}
:class: dropdown
By giving two functions in the plot call, we get the parametric plot $(q(y),y)$ as a function of $y$.
```

```{code-cell}
plot(f, interval..., label="Function",
    legend=:bottomright, grid=:y, xlabel=L"x", ylabel=L"y")
scatter!(x, y, label="initial points")

q = Polynomials.fit(y, x, 2)       # interpolating polynomial
plot!(y -> q(y), y -> y, -0.1, 2.6, l=:dash, label="inverse interpolant")
```

We seek the value of $x$ that makes $y$ zero. This means evaluating $q$ at zero.

```{code-cell}
q(0)
```

Let's restart the process with `BigFloat` numbers to get a convergent sequence.

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
@pt :header=["iteration", "error", "log error", "ratio"] [eachindex(ϵ) ϵ logerr ratios]
```

The convergence is probably superlinear at a rate of $\alpha=1.8$ or so.
