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
[**Demo %s**](#demo-normaleqns-instab)


Because the functions $\sin^2(t)$, $\cos^2(t)$, and $1$ are linearly dependent, we should find that the following matrix is somewhat ill-conditioned.
```{tip}
:class: dropdown
The local variable scoping rule for loops applies to comprehensions as well.
```

```{code-cell}
t = linspace(0, 3, 400)';
A = [ sin(t).^2, cos((1+1e-7)*t).^2, t.^0 ];
kappa = cond(A)
```

Now we set up an artificial linear least-squares problem with a known exact solution that actually makes the residual zero.

```{code-cell}
x = [1; 2; 1];
b = A * x;
```

Using backslash to find the least-squares solution, we get a relative error that is well below $\kappa$ times machine epsilon.

```{code-cell}
x_BS = A \ b;
observed_err = norm(x_BS - x) / norm(x)
max_err = kappa * eps
```

If we formulate and solve via the normal equations, we get a much larger relative error. With $\kappa^2\approx 10^{14}$, we may not be left with more than about 2 accurate digits.

```{code-cell}
N = A'*A;
x_NE = N\(A'*b);
observed_err = norm(x_NE - x) / norm(x)
digits = -log10(observed_err)
```
