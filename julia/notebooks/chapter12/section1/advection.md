---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---

 In the following definition we allow the velocity $c$ to be specified as a parameter in the `ODEProblem`.

```{code-cell}
x, Dₓ, Dₓₓ = FNC.diffper(300, [-4, 4]);
f = (u, c, t) -> -c * (Dₓ*u);
```

The following initial condition isn't mathematically periodic, but the deviation is less than machine precision. We specify RK4 as the solver.  

```{code-cell}
using OrdinaryDiffEq
u_init = @. 1 + exp( -3x^2 )
IVP = ODEProblem(f, u_init, (0., 4.), 2)
sol = solve(IVP, RK4());
```

```{code-cell}
:tags: [hide-input]
using Plots
plt = plot(
    legend=:bottomleft,
    title="Advection with periodic boundary",
    xaxis=(L"x"),  yaxis=(L"u(x,t)"))
for t in (0:4) * 2/3
    plot!(x, sol(t), label=@sprintf("t=%.1f", t))
end
plt
```

An animation shows the solution nicely. The bump moves with speed 2 to the right, reentering on the left as it exits to the right because of the periodic conditions. 

```{code-cell}
:tags: [hide-input]
anim = @animate for t in range(0, 4, 120) 
    plot(x, sol(t),
        title=@sprintf("Advection equation, t = %.2f", t),
        xaxis=(L"x"),  yaxis=([1, 2], L"u(x,t)"),
        dpi=150)
end
mp4(anim, "figures/advection.mp4")
```

