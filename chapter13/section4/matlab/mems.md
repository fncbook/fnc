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
[**Demo %s**](#demo-nonlinear2d-mems)

All we need to define are $\phi$ from {eq}`nonlinpdepde` for the PDE, and a trivial zero function for the boundary condition.

```{code-cell}
lambda = 1.5;
phi = @(X, Y, U, Ux, Uxx, Uy, Uyy) Uxx + Uyy - lambda ./ (U + 1).^2;
g = @(x, y) zeros(size(x));
```

Here is the solution for $m=15$, $n=8$.

```{code-cell}
u = elliptic(phi, g, 15, [0, 2.5], 8, [0, 1]);
disp(sprintf("solution at (2, 0.6) is %.7f", u(2, 0.6)))
```

```{code-cell}
:tags: [hide-input]
x = linspace(0, 2.5, 91);
y = linspace(0, 1, 51);
[mtx, X, Y] = tensorgrid(x, y);
clf,  pcolor(x, y, mtx(u)')
colormap(flipud(sky)),  shading interp,  colorbar
axis equal
xlabel("x"),  ylabel("y")
title("Deflection of a MEMS membrane")
```

In the absence of an exact solution, how can we be confident that the solution is accurate? First, the Levenberg iteration converged without issuing a warning, so we should feel confident that the discrete equations were solved. Assuming that we encoded the PDE correctly, the remaining source of error is truncation from the discretization. We can estimate that by refining the grid a bit and seeing how much the numerical solution changes.

```{code-cell}
x_test = linspace(0, 2.5, 6);
y_test = linspace(0, 1 , 5);
mtx_test = tensorgrid(x_test, y_test);
format long
mtx_test(u)
```

```{code-cell}
u = elliptic(phi, g, 25, [0, 2.5], 14, [0, 1]);
mtx_test(u)
```

The original solution seems to be accurate to about four digits.

