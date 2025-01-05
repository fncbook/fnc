---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---

The equation $u'=(u+t)^2$ gives us some trouble.

```{code-cell}
f(u, p, t) = (t + u)^2
ivp = ODEProblem(f, 1.0, (0.0, 1.0))
sol = solve(ivp, Tsit5());
```

The warning message we received can mean that there is a bug in the formulation of the problem. But if everything has been done correctly, it suggests that the solution may not exist past the indicated time. This is a possibility in nonlinear ODEs.

```{code-cell}
plot(sol, label="";
    xlabel=L"t",  yaxis=(:log10, L"u(t)"),
    title="Finite-time blowup")
```
