---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
include("../../../julia/FNC_init.jl")
```
[**Demo %s**](#demo-wave-speed)

The ODE implementation has to change slightly.

```{code-cell}
ode = function(w,c,t)
    u = extend(w[1:m-1])
    z = w[m:2m]
    du_dt = Dₓ*z
    dz_dt = c.^2 .* (Dₓ*u)
    return [ chop(du_dt); dz_dt ]
end;
```

The variable wave speed is passed as an extra parameter through the IVP solver.

```{code-cell}
using OrdinaryDiffEq
c = @. 1 + (sign(x)+1)/2
IVP = ODEProblem(ode, w_init, (0., 5.), c)
w = solve(IVP, RK4());
```

```{code-cell}
t = range(0, 5, 80)
U = [extend(w(t)[1:m-1]) for t in t]
contour(x, t, hcat(U...)';
    color=:redsblues,  clims=(-1,1),
    levels=24,
    xlabel=L"x",  ylabel=L"t",
    title="Wave equation",
    right_margin=3Plots.mm
    )
```

```{code-cell}
:tags: [hide-input]
anim = @animate for t in range(0,5,181)
    plot(Shape([-1, 0, 0, -1], [-1, -1, 1, 1]), color=RGB(.8, .8, .8), l=0, label="")
    plot!(x, extend(w(t, idxs=1:m-1));
        label=@sprintf("t=%.2f", t), 
        xaxis=(L"x"),  yaxis=([-1, 1], L"u(x,t)"),
        dpi=150,  title="Wave equation, variable speed")
end
mp4(anim, "wave-speed.mp4")
```

Each pass through the interface at $x=0$ generates a reflected and transmitted wave. By conservation of energy, these are both smaller in amplitude than the incoming bump.
