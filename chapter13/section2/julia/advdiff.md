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
[**Demo %s**](#demo-diffadv-advdiff)


The first step is to define a discretization of the domain. 

```{code-cell}
m, n = 50, 36
x, Dx, Dxx = FNC.diffcheb(m, [-1, 1])
y, Dy, Dyy = FNC.diffcheb(n, [-1, 1])
mtx, X, Y, _ = FNC.tensorgrid(x, y)
U₀ = mtx( (x, y) -> (1 + y) * (1 - x)^4 * (1 + x)^2 * (1 - y^4) );
```

There are really two grids now: the full grid and the subset grid of interior points. Since the IVP unknowns are on the interior grid, that is the one we need to change shapes on. We also need the functions `extend` and `chop` to add and remove boundary values.

```{code-cell}
chop = U -> U[2:m, 2:n]
extend = U -> [zeros(m+1) [zeros(1, n-1); U; zeros(1, n-1)] zeros(m+1)]
unvec = u -> reshape(u, m-1, n-1)
pack = U -> vec(chop(U))
unpack = w -> extend(unvec(w))
```

Now we can define and solve the IVP using a stiff solver.

```{code-cell}
using OrdinaryDiffEq
function dw_dt(w, ϵ, t)
    U = unpack(w)
    Ux, Uxx = Dx * U, Dxx * U
    Uyy = U * Dyy'
    du_dt = @. 1 - Ux + ϵ * (Uxx + Uyy)
    return pack(du_dt)
end

IVP = ODEProblem(dw_dt, pack(U₀), (0.0, 2), 0.05)
w = solve(IVP, Rodas4P());
```

When we evaluate the solution at a particular value of $t$, we get a vector of the interior grid values. The same `unpack` function above converts this to a complete matrix of grid values.

```{code-cell}
U = t -> unpack(w(t))
contour(x, y, U(0.5)';
    fill=true,  color=:blues,  levels=20, l=0,
    aspect_ratio=1,  xlabel=L"x",  ylabel=L"y",
    title="Solution at t = 0.5")
```

```{code-cell}
anim = @animate for t in 0:0.02:2
    U = unpack(w(t))
    surface(x, y, U';
        layout=(1, 2),  size=(640, 320),
        xlabel=L"x",  ylabel=L"y",  zaxis=((0, 2), L"u(x,y)"),
        color=:blues,  alpha=0.66,  clims=(0, 2), colorbar=:none,
        title="Advection-diffusion",  dpi=150)
    contour!(x, y, U'; 
        levels=24, 
        aspect_ratio=1,  subplot=2, 
        xlabel=L"x",  ylabel=L"y",
        color=:blues,  clims=(0, 2),  colorbar=:none,
        title=@sprintf("t = %.2f", t))
end
closeall();
mp4(anim, "diffadv-advdiff.mp4");
```

![Advection-diffusion in 2d](diffadv-advdiff.mp4)
