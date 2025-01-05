---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
In {numref}`Demo %s <demo-basics-sing>` we saw an IVP that appears to blow up in a finite amount of time. Because the solution increases so rapidly as it approaches the blowup, adaptive stepping is required even to get close.

```{code-cell}
f(u, p, t) = (t + u)^2
ivp = ODEProblem(f, 1, (0.0, 1.0))
t, u = FNC.rk23(ivp, 1e-5);
```

In fact, the failure of the adaptivity gives a decent idea of when the singularity occurs.

```{code-cell}
plot(t, u;
    legend=:none,
    xlabel=L"t",  yaxis=(:log10, L"u(t)"), 
    title="Finite-time blowup")

tf = t[end]
vline!([tf], l=:dash)
annotate!(tf, 1e5, latexstring(@sprintf("t = %.6f ", tf)), :right)
```
