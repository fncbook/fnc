---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
exec(open("../../../python/FNC_init.py").read())
```
[**Demo %s**](#demo-normaleqns-instab)


Because the functions $\sin^2(t)$, $\cos^2(t)$, and $1$ are linearly dependent, we should find that the following matrix is somewhat ill-conditioned.

```{code-cell}
from numpy.linalg import cond
t = linspace(0, 3, 400)
A = array([ [sin(t)**2, cos((1+1e-7)*t)**2, 1] for t in t ])
kappa = cond(A)
print(f"cond(A) is {kappa:.3e}")
```

Now we set up an artificial linear least-squares problem with a known exact solution that actually makes the residual zero.

```{code-cell}
x = array([1, 2, 1])
b = A @ x
```

Using backslash to find the least-squares solution, we get a relative error that is well below $\kappa$ times machine epsilon.

```{code-cell}
from numpy.linalg import lstsq
x_BS = lstsq(A, b, rcond=None)[0]
print(f"observed error: {norm(x_BS - x) / norm(x):.3e}")
print(f"conditioning bound: {kappa * finfo(float).eps:.3e}")
```

If we formulate and solve via the normal equations, we get a much larger relative error. With $\kappa^2\approx 10^{14}$, we may not be left with more than about 2 accurate digits.

```{code-cell}
N = A.T @ A
x_NE = linalg.solve(N, A.T @ b)
relative_err = norm(x_NE - x) / norm(x)
print(f"observed error: {relative_err:.3e}")
print(f"accurate digits: {-log10(relative_err):.2f}")
```
