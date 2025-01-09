---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
exec(open("../../../python/FNC_init.py").read())
```
[**Demo %s**](#demo-nonlinear2d-allencahn)


The following defines the PDE and a nontrivial Dirichlet boundary condition for the square $[0,1]^2$.

```{code-cell}
phi = lambda x, y, u, ux, uxx, uy, uyy: u * (1 - u**2) + 0.05 * (uxx + uyy)
g = lambda x, y: tanh(5 * (x + 2*y - 1))
```

We solve the PDE and then plot the result.

```{code-cell}
u = FNC.elliptic(phi, g, 36, [0, 1], 36, [0, 1])
```

```{code-cell}
:tags: [hide-input]
x = y = linspace(0, 1, 70)
mtx, X, Y, _, _, _ = FNC.tensorgrid(x, y)
U = mtx(u)
pcolormesh(X.T, Y.T, U.T, cmap="viridis")
xlabel("$x$"),  ylabel("$y$"),  axis("equal")
colorbar(),  title("Steady Allen–Cahn equation");
```
