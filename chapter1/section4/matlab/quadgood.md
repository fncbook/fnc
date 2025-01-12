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
[**Demo %s**](#demo-stability-quadgood)

We repeat the rootfinding experiment of {numref}`Demo %s <demo-stability-quadbad>` with an alternative algorithm.

```{code-cell}
a = 1;  b = -(1e6 + 1e-6);  c = 1;
```

First, we find the "good" root using the quadratic formula.

```{code-cell}
x1 = (-b + sqrt(b^2 - 4*a*c)) / (2*a);
```

Then we use the identity $x_1x_2=\frac{c}{a}$ to compute the smaller root:

```{code-cell}
x2 = c / (a*x1)
```

This matches the exact root to the displayed digits; to be sure we have an accurate result, we compute its relative error.

```{code-cell}
abs(x2 - 1e-6) / 1e-6
```
