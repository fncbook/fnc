---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
We study the convergence of AB4 using the IVP $u'=\sin[(u+t)^2]$ over $0\le t \le 4$, with $u(0)=-1$. As usual, `solve` is called to give an accurate reference solution.

```{code-cell}
using OrdinaryDiffEq
ivp = ODEProblem((u, p, t) -> sin((t + u)^2), -1.0, (0.0, 4.0))
u_ref = solve(ivp, Tsit5(), reltol=1e-14, abstol=1e-14);
```

Now we perform a convergence study of the AB4 code.

```{code-cell}
n = @. [round(Int, 4 * 10^k) for k in 0:0.5:3]
err = []
for n in n
    t, u = FNC.ab4(ivp, n)
    push!(err, norm(u_ref.(t) - u, Inf))
end
@pt :header=["n", "inf-norm error"] [n err]
```

The method should converge as $O(h^4)$, so a log-log scale is appropriate for the errors.

```{code-cell}
using Plots
plot(n, err, m=3, 
    label="AB4",  legend=:bottomleft,
    xaxis=(:log10, L"n"),  yaxis=(:log10, "inf-norm error"),
    title="Convergence of AB4")

plot!(n, 0.1 * err[end] * (n / n[end]) .^ (-4), l=:dash, label=L"O(n^{-4})")
```
