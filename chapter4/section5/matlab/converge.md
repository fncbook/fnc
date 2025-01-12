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
[**Demo %s**](#demo-newtonsys-converge)

A system of nonlinear equations is defined by its residual and Jacobian.
```{tip}
:class: dropdown
This function needs to be defined within a script file or in a file of its own with the `.m` extension.
```

```{literalinclude} f45_nlsystem.m
:language: matlab
```

Since our system function is defined in an external file here, we need to use `@` in order to reference it as a function argument. 

```{code-cell}
nlsystem = @f45_nlsystem;
x1 = [0; 0; 0];    % column vector!
x = newtonsys(nlsystem, x1);
num_iter = size(x, 2)
```

Let's compute the residual of the last result in order to check the quality.

```{code-cell}
r = x(:, end)
back_err = norm(nlsystem(r))
```

We take the sequence errors in the first component of the solution, applying the log so that we can look at the exponents.

```{code-cell}
log10( abs(x(1, 1:end-1) - r(1)) )'
```

This sequence looks to be nearly doubling at each iteration, which is a good sign of quadratic convergence.
