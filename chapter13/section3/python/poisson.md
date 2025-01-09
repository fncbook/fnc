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
[**Demo %s**](#demo-laplace-poisson)


First we define the problem on $[0,1]\times[0,2]$.

```{code-cell}
f = lambda x, y: -sin(3 * x * y - 4 * y) * (9 * y**2 + (3 * x - 4) ** 2)
g = lambda x, y: sin(3 * x * y - 4 * y)
xspan = [0, 1]
yspan = [0, 2]
```

Here is the finite-difference solution.

```{code-cell}
U, X, Y = FNC.poissonfd(f, g, 50, xspan, 80, yspan)
```

```{code-cell}
:tags: [hide-input]
pcolormesh(X.T, Y.T, U.T, cmap="viridis")
xlabel("$x$"),  ylabel("$y$"),  axis("equal")
colorbar(), title("Solution of Poisson's equation");
```

Since this is an artificial problem with a known solution, we can plot the error, which is a smooth function of $x$ and $y$. It must be zero on the boundary; otherwise, we have implemented boundary conditions incorrectly.

```{code-cell}
error = g(X, Y) - U    # because we set up g as the exact solution
M = max(abs(error))

pcolormesh(X.T, Y.T, error.T, vmin=-M, vmax=M, cmap="RdBu")
xlabel("$x$"),  ylabel("$y$"),  axis("equal")
colorbar(),  title("Error");
```
