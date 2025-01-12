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
[**Demo %s**](#demo-nonlinear-advdiff)


```{code-cell}
phi = @(X, Y, U, Ux, Uxx, Uy, Uyy) 1 - Ux - 2*Uy + 0.05*(Uxx + Uyy);
g = @(x, y) zeros(size(x));
u = elliptic(phi, g, 32, [-1, 1], 32, [-1, 1]);
```

```{code-cell}
:tags: [hide-input]
x = linspace(-1, 1, 80);
y = x;
mtx = tensorgrid(x, y);
clf,  pcolor(x, y, mtx(u)')
colormap(parula),  shading interp,  colorbar
axis equal,  xlabel("x"),  ylabel("y")
title("Steady advection–diffusion")
```

