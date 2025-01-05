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
[**Demo %s**](#demo-secant-converge)

We check the convergence of the secant method from {numref}`Demo %s <demo-secant-line>`. Again we will use extended precision to get a longer sequence than double precision allows.

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
@pt :header=["iteration", "error", "log error", "ratio"] [eachindex(ϵ) ϵ logerr ratios]
```

As expected, this settles in at around 1.618.
