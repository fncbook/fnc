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
[**Demo %s**](#demo-nonlinear-mems)


Here is the problem definition. We use a truncated domain to avoid division by zero at $r=0$.

```{code-cell}
lambda = 0.5;
phi = @(r,w,dwdr) lambda./w.^2 - dwdr./r;
ga = @(w, dw) dw;
gb = @(w, dw) w - 1;
a = eps;  b = 1;
```

First we try a constant function as the initialization.

```{code-cell}
init = ones(301, 1);
[r, w1] = bvp(phi, a, b, ga, gb, init);

clf,  plot(r, w1)
xlabel('r'),  ylabel('w(r)')
title('Solution of the MEMS BVP')
```

It's not necessary that the initialization satisfy the boundary conditions. In fact, by choosing a different constant function as the initial guess, we arrive at another valid solution.

```{code-cell}
init = 0.5 * ones(301, 1);
[r, w2] = bvp(phi, a, b, ga, gb, init);
hold on,  plot(r, w2)
title("Two solutions of the MEMS BVP")
```
