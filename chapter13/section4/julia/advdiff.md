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
[**Demo %s**](#demo-nonlinear-advdiff)


```{code-cell}
ϕ = (X, Y, U, Ux, Uxx, Uy, Uyy) -> @. 1 - Ux - 2Uy + 0.05 * (Uxx + Uyy)
g = (x, y) -> 0
u = FNC.elliptic(ϕ, g, 32, [-1, 1], 32, [-1, 1]);
```

```{code-cell}
x = y = range(-1, 1, 80)
U = [u(x, y) for x in x, y in y]
contourf(x, y, U';
    color=:viridis, 
    aspect_ratio=1,
    xlabel=L"x",  ylabel=L"y",  zlabel=L"u(x,y)",
    title="Steady advection–diffusion")
```

