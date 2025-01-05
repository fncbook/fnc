---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---

```{code-cell}
K = 3;  σ = 0.06;  r = 0.08;  Smax = 8;
ϕ = (t, x, u, uₓ, uₓₓ) -> σ^2/2 * (x^2 * uₓₓ) + r*x*uₓ - r*u
g₁ = (u, uₓ) -> u
g₂ = (u, uₓ) -> uₓ - 1;
```

```{code-cell}
u₀ = x -> max(0, x - K)
x, u = FNC.parabolic(ϕ, (0, Smax), 80, g₁, g₂, (0, 15), u₀);
```

```{code-cell}
:tags: [hide-input]
anim = @animate for t in range(0, 15, 151) 
    plot(x, u(t);
        label=@sprintf("t=%.4f", t),  legend=:topleft,
        xaxis=(L"x"),  yaxis=(L"u(x,t)", (-0.5, 8)), 
        dpi=150,  title="Black–Scholes equation")
end
mp4(anim, "figures/boundaries-bs.mp4", fps=30)
```

Recall that $u$ is the value of the call option, and time runs backward from the strike time. The longer the horizon, the more value the option has due to anticipated growth in the stock price.
