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
[**Demo %s**](#demo-adapt-basic)

Let's run adaptive RK on  $u'=e^{t-u\sin u}$.

```{code-cell}
using OrdinaryDiffEq, Plots
f(u, p, t) = exp(t - u * sin(u))
ivp = ODEProblem(f, 0, (0.0, 5.0))
t, u = FNC.rk23(ivp, 1e-5)
plot(t, u, m=2,
    xlabel=L"t",  ylabel=L"u(t)", 
    title="Adaptive IVP solution")
```

The solution makes a very abrupt change near $t=2.4$. The resulting time steps vary over three orders of magnitude.

```{code-cell}
Δt = diff(t)
plot(t[1:end-1], Δt;
    xaxis=(L"t", (0, 5)), yaxis=(:log10, "step size"),
    title="Adaptive step sizes")
```

If we had to run with a uniform step size to get this accuracy, it would be

```{code-cell}
println("minimum step size = $(minimum(Δt))")
```

On the other hand, the average step size that was actually taken was

```{code-cell}
println("average step size = $(sum(Δt)/(length(t)-1))")
```

We took fewer steps by a factor of almost 1000! Even accounting for the extra stage per step and the occasional rejected step, the savings are clear.

