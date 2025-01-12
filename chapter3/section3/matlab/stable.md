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
[**Demo %s**](#demo-qr-stable)

We'll repeat the experiment of {numref}`Demo {number} <demo-normaleqns-instab>`, which exposed instability in the normal equations. 

```{code-cell}
t = linspace(0, 3, 400)';
A = [ sin(t).^2, cos((1+1e-7)*t).^2, t.^0 ];
x = [1; 2; 1];
b = A * x;
```

The error in the solution by {numref}`Function {number} <function-lsqrfact>` is similar to the bound predicted by the condition number.

```{code-cell}
observed_error = norm(lsqrfact(A, b) - x) / norm(x)
error_bound = cond(A) * eps
```
