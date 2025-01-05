---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
include("../../../julia/FNC_init.jl")
```
[**Demo %s**](#demo-stability-rungefix)

Here again is the function from {numref}`Demo {number} <demo-stability-runge>` that provoked the Runge phenomenon when using equispaced nodes.

```{code-cell} 
f(x) = 1 / (x^2 + 16);
```

```{code-cell} 
:tags: [hide-input]

plot(label="", xaxis=(L"x"), yaxis=(:log10, L"|f(x)-p(x)|", [1e-20, 1]))
x = range(-1, 1, 2001)
for (k, n) in enumerate([4, 10, 16, 40])
    t = [-cospi(k / n) for k in 0:n]
    y = f.(t)                           # interpolation data
    p = FNC.polyinterp(t, y)
    err = @. abs(f(x) - p(x))
    plot!(x, err, m=(1, :o, stroke(0)), label="degree $n")
end
title!("Error for Chebyshev interpolants")
```

By degree 16 the error is uniformly within machine epsilon, and, importantly, it stays there as $n$ increases. Note that as predicted by the error indicator function, the error is uniform over the interval at each value of $n$.
