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
[**Demo %s**](#demo-nonlinear-advdiff)


```{code-cell}
phi = lambda x, y, u, ux, uxx, uy, uyy: 1 - ux - 2*uy + 0.05 * (uxx + uyy)
g = lambda x, y: 0
u = FNC.elliptic(phi, g, 32, [-1, 1], 32, [-1, 1])
```

```{code-cell}
x = y = linspace(-1, 1, 70)
mtx, X, Y, _, _, _ = FNC.tensorgrid(x, y)
U = mtx(u)

pcolormesh(X.T, Y.T, U.T, cmap="viridis")
xlabel("$x$"),  ylabel("$y$"),  axis("equal")
colorbar()
title("Steady advection–diffusion");
```

