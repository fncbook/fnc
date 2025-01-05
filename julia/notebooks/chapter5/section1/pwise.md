---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering: 
  headings: false
---
Let us recall the data from {numref}`Demo %s <demo-interpolation-global>`.

```{code-cell}
n = 12
t = range(-1, 1, n+1)
y = @. t^2 + t + 0.5 * sin(20t)
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
using Dierckx
p = Spline1D(t, y)
plot!(x -> p(x), -1, 1, label="piecewise cubic")
```
