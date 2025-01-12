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
[**Demo %s**](#demo-nonlinear2d-allencahn)


The following defines the PDE and a nontrivial Dirichlet boundary condition for the square $[0,1]^2$.

```{code-cell}
phi = @(X, Y, U, Ux, Uxx, Uy, Uyy) U .* (1-U.^2) + 0.05*(Uxx + Uyy);
g = @(x, y) tanh(5*(x + 2*y - 1));
```

We solve the PDE and then plot the result.

```{code-cell}
u = elliptic(phi, g, 36, [0, 1], 36, [0, 1]);
```

```{code-cell}
:tags: [hide-input]
x = linspace(0, 1, 80);
y = x;
mtx = tensorgrid(x, y);
clf,  pcolor(x, y, mtx(u)')
colormap(parula),  shading interp,  colorbar
axis equal,  xlabel("x"),  ylabel("y")
title("Steady Allen–Cahn")
```
