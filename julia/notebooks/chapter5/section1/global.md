---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering: 
  headings: false
---
Here are some points that we could consider to be observations of an unknown function on $[-1,1]$.

```{code-cell}
using Plots
n = 5
t = range(-1, 1, n+1)
y = @. t^2 + t + 0.05 * sin(20t)
scatter(t, y, label="data", legend=:top)
```

```{index} ! Julia; fit
```

The polynomial interpolant, as computed using `fit`, looks very sensible. It's the kind of function you'd take home to meet your parents.

```{code-cell}
using Polynomials
p = Polynomials.fit(t, y, n)     # interpolating polynomial
plot!(p, -1, 1, label="interpolant")
```

But now consider a different set of points generated in almost exactly the same way.

```{code-cell}
n = 18
t = range(-1, 1, n+1)
y = @. t^2 + t + 0.05 * sin(20t)
scatter(t, y, label="data", leg=:top)
```

The points themselves are unremarkable. But take a look at what happens to the polynomial interpolant.

```{code-cell}
p = Polynomials.fit(t, y, n)
x = range(-1, 1, 1000)    # use a lot of points
plot!(x, p.(x), label="interpolant")
```

Surely there must be functions that are more intuitively representative of those points!
