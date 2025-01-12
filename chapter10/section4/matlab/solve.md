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
[**Demo %s**](#demo-linear-solve)


```{code-cell}
exact = @(x) exp(sin(x));
```

The problem is presented above in our standard form, so we can identify the coefficient functions in the ODE. Each should be coded as a function.

```{code-cell}
p = @(x) -cos(x);
q = @(x) sin(x);
r = @(x) 0*x;      % not a scalar 
```

We solve the BVP and compare the result to the exact solution.

```{code-cell}
[x, u] = bvplin(p, q, r, 0, 3*pi/2, 1, exp(-1), 30);
```

```{code-cell}
clf,  subplot(2, 1, 1)
plot(x, u)
ylabel('solution')
title('Solution of a linear BVP')
subplot(2, 1, 2)
plot(x, exact(x) - u, 'o-')
ylabel('error')
```
