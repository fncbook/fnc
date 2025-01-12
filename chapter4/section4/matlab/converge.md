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
[**Demo %s**](#demo-secant-converge)

We check the convergence of the secant method from {numref}`Demo %s <demo-secant-line>`. 

```{code-cell}
f = @(x) x .* exp(x) - 2;
x = secant(f, 1, 0.5);
```

We don't know the exact root, so we use `fzero` to get a good proxy.

```{code-cell}
r = fzero(f, 1);
```

Here is the sequence of errors.

```{code-cell}
format short e
err = r - x(1:end-1)'
```

It's not easy to see the convergence rate by staring at these numbers. We can use {eq}`superlinear-rate` to try to expose the superlinear convergence rate.

```{code-cell}
logerr = log(abs(err));
ratios = logerr(2:end) ./ logerr(1:end-1)
```

As expected, this settles in at around 1.618.
