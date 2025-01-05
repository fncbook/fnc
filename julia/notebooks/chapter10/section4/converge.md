---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{code-cell}
λ = 10
exact = x -> sinh(λ * x) / sinh(λ) - 1;
```

The following functions define the ODE.

```{code-cell}
p = x -> 0
q = x -> -λ^2
r = x -> λ^2;
```

We compare the computed solution to the exact one for increasing $n$.

```{code-cell}
n = 5 * [round(Int, 10^d) for d in 0:0.25:3]
err = zeros(length(n))
for (k, n) in enumerate(n)
    x, u = FNC.bvplin(p, q, r, [0, 1], -1, 0, n)
    err[k] = norm(exact.(x) - u, Inf)
end
data = (n = n[1:4:end], err = err[1:4:end])
@pt :header = ["n", "inf-norm error"] data
```

Each factor of 10 in $n$ reduces error by a factor of 100, which is indicative of second-order convergence.

```{code-cell}
plot(n, err, m = :o,
    label = "observed",
    xaxis = (:log10, L"n"),
    yaxis = (:log10, "inf-norm error"),
    title = "Convergence for a linear BVP")
plot!(n, 0.25 * n .^ (-2), l = (:dash, :gray), label = "2nd order")
```
