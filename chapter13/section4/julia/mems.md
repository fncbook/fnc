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
[**Demo %s**](#demo-nonlinear2d-mems)

All we need to define are $\phi$ from {eq}`nonlinpdepde` for the PDE, and a trivial zero function for the boundary condition.

```{code-cell}
λ = 1.5
ϕ = (X, Y, U, Ux, Uxx, Uy, Uyy) -> @. Uxx + Uyy - λ / (U + 1)^2;
g = (x, y) -> 0;
```

Here is the solution for $m=15$, $n=8$.

```{code-cell}
u = FNC.elliptic(ϕ, g, 15, [0, 2.5], 8, [0, 1]);
```

```{code-cell}
x = range(0, 2.5, 100)
y = range(0, 1, 50)
U = [u(x, y) for x in x, y in y]
contourf(x, y, U';
    color=:blues,  l=0,
    aspect_ratio=1,
    xlabel=L"x",  ylabel=L"y",  zlabel=L"u(x,y)",
    title="Deflection of a MEMS membrane",
    right_margin=3Plots.mm)
```

In the absence of an exact solution, how can we be confident that the solution is accurate? First, the Levenberg iteration converged without issuing a warning, so we should feel confident that the discrete equations were solved. We can check the boundary values easily. For example,

```{code-cell}
x_test = range(0, 2.5, 100)
norm([u(x, 0) - g(x, 0) for x in x_test], Inf)
```

Assuming that we encoded the PDE correctly, the remaining source error is truncation from the discretization. We can estimate that by refining the grid a bit and seeing how much the numerical solution changes.

```{code-cell}
x_test = range(0, 2.5, 6)
y_test = range(0, 1, 6)
mtx_test, _ = FNC.tensorgrid(x_test, y_test)
mtx_test(u)
```

```{code-cell}
u = FNC.elliptic(ϕ, g, 25, [0, 2.5], 14, [0, 1]);
mtx_test(u)
```

The original solution seems to be accurate to about four digits.

