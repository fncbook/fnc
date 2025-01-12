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
[**Demo %s**](#demo-newton-converge)

We again look at finding a solution of $x e^x=2$ near $x=1$. To apply Newton's method, we need to calculate values of both the residual function $f$ and its derivative.

```{code-cell}
f = @(x) x.*exp(x) - 2;
dfdx = @(x) exp(x).*(x+1);
```

We don't know the exact root, so we use `nlsolve` to determine a proxy for it.

```{code-cell}
format long,  r = fzero(f,1)
```

We use $x_1=1$ as a starting guess and apply the iteration in a loop, storing the sequence of iterates in a vector.

```{code-cell}
x = 1;
for k = 1:6
    x(k+1) = x(k) - f(x(k)) / dfdx(x(k));
end
x
```

Here is the sequence of errors.

```{code-cell}
format short e
err = x' - r
```

The exponents in the scientific notation definitely suggest a squaring sequence. We can check the evolution of the ratio in {eq}`quadratictest`.

```{code-cell}
format short
logerr = log(abs(err))
```

The clear convergence to 2 above constitutes good evidence of quadratic convergence.
