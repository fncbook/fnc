---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{code-cell}
f(t) = π * sqrt( cospi(t)^2 + sinpi(t)^2 / 4 );
n = 4:4:48
perim = zeros(size(n))
for (k, n) in enumerate(n)
    h = 2 / n
    t = @. h * (0:n-1) - 1
    perim[k] = h * sum(f.(t))
end
err = @. abs(perim - perim[end])    # use last value as "exact"
@ptconf formatters=ft_printf(["%d", "%.15f", "%.2e"], 1:3)
@pt :header=["n", "perimeter", "error estimate"] [n perim err][1:end-1, :]
```
The approximations gain about one digit of accuracy for each constant increment of $n$, which is consistent with spectral convergence.
