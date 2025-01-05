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
[**Demo %s**](#demo-stability-errcheb)

Now we look at the error indicator function $\Phi$ for Chebyshev node sets.

```{code-cell} 
:tags: [hide-input]
plot(xaxis=(L"x"), yaxis=(:log10, L"|\Phi(x)|", [1e-18, 1e-2]))
x = range(-1, 1, 2001)
for n in 10:10:50
    t = [-cospi(k / n) for k in 0:n]
    Φ(x) = prod(x - t for t in t)
    plot!(x, abs.(Φ.(x)), m=(1, :o, stroke(0)), label="n=$n")
end
title!("Error indicator for Chebyshev nodes")
```

In contrast to the equispaced case, $|\Phi|$ decreases exponentially with $n$ almost uniformly across the interval.
