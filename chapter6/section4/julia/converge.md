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
[**Demo %s**](#demo-rk-converge)

We solve the IVP $u'=\sin[(u+t)^2]$ over $0\le t \le 4$, with $u(0)=-1$.

```{code-cell}
using OrdinaryDiffEq
f(u, p, t) = sin((t + u)^2)
tspan = (0.0, 4.0)
u₀ = -1.0
ivp = ODEProblem(f, u₀, tspan)
```

We use a `DifferentialEquations` solver to construct an accurate approximation to the exact solution.

```{code-cell}
u_ref = solve(ivp, Tsit5(), reltol=1e-14, abstol=1e-14);
```

Now we perform a convergence study of our two Runge–Kutta implementations.

```{code-cell}
n = [round(Int, 2 * 10^k) for k in 0:0.5:3]
err = zeros(length(n), 2)
for (k, n) in enumerate(n)
    t, u = FNC.ie2(ivp, n)
    err[k, 1] = norm(u_ref.(t) - u, Inf)
    t, u = FNC.rk4(ivp, n)
    err[k, 2] = norm(u_ref.(t) - u, Inf)
end
@pt :header=["n", "IE2 error", "RK4 error"] [n err]
```

The amount of computational work at each time step is assumed to be proportional to the number of stages. Let's compare on an apples-to-apples basis by using the number of $f$-evaluations on the horizontal axis.

```{code-cell}
using Plots
plot([2n 4n], err;
    m=3, label=["IE2" "RK4"], legend=:bottomleft,
    xaxis=(:log10, "f-evaluations"),  yaxis=(:log10, "inf-norm error"),
    title="Convergence of RK methods")

plot!(2n, 0.1 * err[end,1] * (n / n[end]) .^ (-2), l=:dash, label=L"O(n^{-2})")
plot!(4n, 0.1 * err[end,2] * (n / n[end]) .^ (-4), l=:dash, label=L"O(n^{-4})")
```

The fourth-order variant is more efficient in this problem over a wide range of accuracy.
