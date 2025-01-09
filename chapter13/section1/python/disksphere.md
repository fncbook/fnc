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
[**Demo %s**](#demo-tensorprod-disksphere)

For a function given in polar form, such as $f(r,\theta)=1-r^4$, construction of a function over the unit disk is straightforward using a grid in $(r,\theta)$ space.

```{code-cell}
r = linspace(0, 1, 41)
theta = linspace(0, 2*pi, 121)
mtx, R, Theta, _, _, _ = FNC.tensorgrid(r, theta)

F = mtx(lambda r, theta: 1 - r**4)    

contourf(R.T, Theta.T, F.T, levels=20)
colorbar()
xlabel("$r$"),  ylabel("$\\theta$");
```

Of course, we are used to seeing such plots over the $(x,y)$ plane, not the $(r,\theta)$ plane. For this we create matrices for the coordinate functions $x$ and $y$.

```{code-cell}
X, Y = R * cos(Theta), R * sin(Theta)
contourf(X.T, Y.T, F.T, levels=20)
colorbar(),  axis("equal")
xlabel("$x$"),  ylabel("$y$");
```

In such functions the values along the line $r=0$ must be identical, and the values on the line $\theta=0$ should be identical to those on $\theta=2\pi$. Otherwise the interpretation of the domain as the unit disk is nonsensical. If the function is defined in terms of $x$ and $y$, then those can be defined in terms of $r$ and $\theta$ using {eq}`unitdiskparam`.

