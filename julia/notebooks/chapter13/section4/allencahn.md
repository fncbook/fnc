---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---

The following defines the PDE and a nontrivial Dirichlet boundary condition for the square $[0,1]^2$.

```{code-cell}
ϕ = (X, Y, U, Ux, Uxx, Uy, Uyy) -> @. U * (1 - U^2) + 0.05 * (Uxx + Uyy)
g = (x, y) -> tanh(5 * (x + 2y - 1));
```

We solve the PDE and then plot the result.

```{code-cell}
u = FNC.elliptic(ϕ, g, 36, [0, 1], 36, [0, 1]);
```

```{code-cell}
x = y = range(0, 1, 80)
U = [u(x, y) for x in x, y in y]
contourf(x, y, U';
    color=:viridis, 
    aspect_ratio=1,
    xlabel=L"x",  ylabel=L"y",  zlabel=L"u(x,y)", 
    title="Steady Allen-Cahn",
    right_margin=3Plots.mm)
```
