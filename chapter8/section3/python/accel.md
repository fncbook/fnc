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
[**Demo %s**](#demo-inviter-accel)

```{code-cell}
ev = array([1, -0.75, 0.6, -0.4, 0])
A = triu(ones([5, 5]), 1) + diag(ev)    # triangular matrix, eigs on diagonal
```

We begin with a shift $s=0.7$, which is closest to the eigenvalue 0.6.

```{code-cell}
from numpy.linalg import solve
s = 0.7
x = ones(5)
y = solve(A - s * eye(5), x)
beta = x[0] / y[0] + s
print(f"latest estimate: {beta:.8f}")
```

Note that the result is not yet any closer to the targeted 0.6. But we proceed (without being too picky about normalization here).

```{code-cell}
s = beta
x = y / y[0]
y = solve(A - s * eye(5), x)
beta = x[0] / y[0] + s
print(f"latest estimate: {beta:.8f}")
```

Still not much apparent progress. However, in just a few more iterations the results are dramatically better.

```{code-cell}
for k in range(4):
    s = beta
    x = y / y[0]
    y = solve(A - s * eye(5), x)
    beta = x[0] / y[0] + s
    print(f"latest estimate: {beta:.12f}")
```
