---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
Here is a $4\times 4$ Chebyshev differentiation matrix.

```{code-cell}
t, Dₓ = FNC.diffcheb(3, [-1, 1])
Dₓ
```

We again test the convergence rate.

```{code-cell}
f = x -> x + exp(sin(4 * x));
df_dx = x -> 1 + 4 * exp(sin(4 * x)) * cos(4 * x);
d2f_dx2 = x -> 4 * exp(sin(4 * x)) * (4 * cos(4 * x)^2 - 4 * sin(4 * x));
```

```{code-cell}
n = 5:5:70
err1 = zeros(size(n))
err2 = zeros(size(n))
for (k, n) in enumerate(n)
    t, Dₓ, Dₓₓ = FNC.diffcheb(n, [-1, 1])
    y = f.(t)
    err1[k] = norm(df_dx.(t) - Dₓ * y, Inf)
    err2[k] = norm(d2f_dx2.(t) - Dₓₓ * y, Inf)
end
```

Since we expect a spectral convergence rate, we use a semi-log plot for the error.

```{code-cell}
plot(n, [err1 err2]; m = :o,
    label=[L"f'" L"f''"],
    xaxis=(L"n"),  yaxis = (:log10, "max error"),
    title="Convergence of Chebyshev derivatives")
```
