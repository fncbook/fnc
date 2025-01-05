---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{code-cell}
:tags: [hide-input]
using Plots
f(x) = 1 / (1 + x^2)
plot(f, -4, 4, layout=(2, 1),
    xlabel=L"x", 
    yaxis=(:log10, L"f(x)", (1e-16, 2)),
    title="Original integrand")

ξ(t) = sinh( π * sinh(t) / 2 )
dξ_dt(t) = π/2 * cosh(t) * cosh(π * sinh(t) / 2)
g(t) = f(ξ(t)) * dξ_dt(t)

plot!(g,-4, 4, subplot=2,
    xlabel=L"t",
    yaxis=(:log10, L"f(x(t))\cdot x'(t)", (1e-16, 2)),
    title="Transformed integrand")
```

This graph suggests that we capture all of the integrand values that are larger than machine epsilon by integrating in $t$ from $-4$ to $4$.
