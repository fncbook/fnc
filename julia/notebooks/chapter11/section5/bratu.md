---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---

```{code-cell}
ϕ = (t, x, u, uₓ, uₓₓ) -> u^2 + uₓₓ
g₁ = (u, uₓ) -> u
g₂ = (u, uₓ) -> uₓ
init = x -> 400x^4 * (1 - x)^2
x, u = FNC.parabolic(ϕ, (0, 1), 60, g₁, g₂, (0, 0.1), init);
```

```{code-cell}
:tags: [hide-input]
anim = @animate for t in range(0, 0.1, length=101) 
    plot(x, u(t);
        label=@sprintf("t=%.4f", t),  legend=:topleft,
        xaxis=(L"x"),  yaxis=(L"u(x,t)", (0, 10)),
        dpi=150, title="Heat equation with source")
end
mp4(anim, "figures/boundaries-source.mp4", fps=30)
```
