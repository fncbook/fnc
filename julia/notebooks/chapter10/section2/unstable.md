---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---

```{code-cell}
plt = plot(
    xaxis = (L"x"),
    yaxis = ([-1.2, 0.5], L"u(x)"),
    title = "Shooting instability",
    leg = :topleft,
)
for λ in 6:4:18
    g₁(u, du) = u + 1
    g₂(u, du) = u
    ϕ = (x, u, du_dx) -> λ^2 * (u + 1)
    x, u = FNC.shoot(ϕ, (0.0, 1.0), g₁, g₂, [-1, 0])
    plot!(x, u, label = "λ=$λ")
end
plt
```

The numerical solutions evidently don't satisfy the right boundary condition as $\lambda$ increases, which makes them invalid. 
