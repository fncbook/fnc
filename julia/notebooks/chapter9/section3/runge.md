---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
This function has infinitely many continuous derivatives on the entire real line and looks easy to approximate over $[-1,1]$.

```{code-cell} 
f(x) = 1 / (x^2 + 16)
plot(f, -1, 1, title="Test function", legend=:none)
```

We start by doing equispaced polynomial interpolation for some small values of $n$.

```{code-cell} 
:tags: [hide-input]
plot(xaxis=(L"x"), yaxis=(:log10, L"|f(x)-p(x)|", [1e-20, 1]))
x = range(-1, 1, 2501)
n = 4:4:12
for (k, n) in enumerate(n)
    t = range(-1, 1, n+1)           # equally spaced nodes
    y = f.(t)                       # interpolation data
    p = FNC.polyinterp(t, y)
    err = @. abs(f(x) - p(x))
    plot!(x, err, m=(1, :o, stroke(0)), label="degree $n")
end
title!("Error for low degrees")
```

The convergence so far appears rather good, though not uniformly so. However, notice what happens as we continue to increase the degree.

```{code-cell} 
:tags: [hide-input]
n = @. 12 + 15 * (1:3)
plot(xaxis=(L"x"), yaxis=(:log10, L"|f(x)-p(x)|", [1e-20, 1]))
for (k, n) in enumerate(n)
    t = range(-1, 1, n+1)           # equally spaced nodes
    y = f.(t)                       # interpolation data
    p = FNC.polyinterp(t, y)
    err = @. abs(f(x) - p(x))
    plot!(x, err, m=(1, :o, stroke(0)), label="degree $n")
end
title!("Error for higher degrees")
```

The convergence in the middle can't get any better than machine precision relative to the function values. So maintaining the growing gap between the center and the ends pushes the error curves upward exponentially fast at the ends, wrecking the convergence.
