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
[**Demo %s**](#demo-galerkin-fem)


Here are the coefficient function definitions. Even though $s$ is a constant, it has to be defined as a function for {numref}`Function {number} <function-fem>` to use it.

```{code-cell}
c = @(x) x.^2;
q = @(x) 4 * ones(size(x));
f = @(x) sin(pi * x);
```

```{code-cell}
[x,u] = fem(c, q, f, 0, 1, 50);
clf,  plot(x, u)
xlabel('x'),  ylabel('u')
title('Solution by finite elements')
```
