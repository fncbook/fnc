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
[**Demo %s**](#demo-quasi-levenberg)


To solve a nonlinear system, we need to code only the function defining the system, and not its Jacobian.
```{tip}
:class: dropdown
A rule of thumb is that if you use a function as an input argument for another function, there needs to be an `@` involved once: either for an anonymous definition or to reference a function defined elsewhere. 
```{literalinclude} f45_nlsystem.m
:language: matlab
```
In all other respects usage is the same as for the `newtonsys` function.
```{code-cell}
f = @f46_nlsystem;
x1 = [0; 0; 0];   
x = levenberg(f, x1);
```
It's always a good idea to check the accuracy of the root, by measuring the residual (backward error).
```{code-cell}
r = x(:, end)
backward_err = norm(f(r))
```
Looking at the convergence of the first component, we find a rate between linear and quadratic, like with the secant method.
```{code-cell}
log10( abs(x(1, 1:end-1) - r(1)) )'
```
