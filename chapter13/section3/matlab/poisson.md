---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
cd  /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init
```
[**Demo %s**](#demo-laplace-poisson)


First we define the problem on $[0,1]\times[0,2]$.

```{code-cell}
f = @(x, y) -sin(3*x .* y - 4*y) .* (9*y.^2 + (3*x - 4).^2);
g = @(x, y) sin(3*x .* y - 4*y);
xspan = [0, 1];
yspan = [0, 2];
```

Here is the finite-difference solution.

```{code-cell}
[X, Y, U] = poissonfd(f, g, 40, xspan, 60, yspan);

clf, surf(X', Y', U')
colormap(parula),  shading interp
colorbar
title("Solution of Poisson's equation")
xlabel("x"),  ylabel("y"),  zlabel("u(x,y)")
```

Since this is an artificial problem with a known solution, we can plot the error, which is a smooth function of $x$ and $y$. It must be zero on the boundary; otherwise, we have implemented boundary conditions incorrectly.

```{code-cell}
err = g(X, Y) - U;
mx = max(abs(vec(err)));
pcolor(X', Y', err')
colormap(redsblues),  shading interp
clim([-mx, mx]),  colorbar
axis equal,  xlabel("x"),  ylabel("y")
title("Error")
```
