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
[**Demo %s**](#demo-tensorprod-gridfun)

Here is the grid from {numref}`Example {number} <example-tensorprod-smallgrid>`.

```{code-cell}
m = 4
x = linspace(0, 2, m+1)
n = 2
y = linspace(1, 3, n+1)
```

For a given $f(x,y)$ we can find $\operatorname{mtx}(f)$ by using a comprehension syntax.

```{code-cell}
f = lambda x, y: cos(pi * x * y - y)
F = array( [ [f(xi, yj) for yj in y] for xi in x ] )
print(F)
```

We can make a nice plot of the function by first choosing a much finer grid. However, the contour and surface plotting functions expect the *transpose* of mtx($f$).
```{tip}
:class: dropdown
To emphasize departures from a zero level, use a colormap such as `RdBu` and set the color limits to be symmetric around zero.
```

::::{warning}
The contour and surface plotting functions expect the *transpose* of the outputs of `mtx`. If you forget to do that, the $x$ and $y$ axes will be swapped.
::::

```{code-cell}
m, n = 80, 70
x = linspace(0, 2, m+1)
y = linspace(1, 3, n+1)
mtx, X, Y, _, _, _ = FNC.tensorgrid(x, y)
F = mtx(f)

pcolormesh(X.T, Y.T, F.T, cmap="RdBu", vmin=-1, vmax=1, shading="gouraud")
axis("equal"),  colorbar()
xlabel("$x$"),  ylabel("$y$");
```
