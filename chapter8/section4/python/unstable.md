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
[**Demo %s**](#demo-subspace-unstable)

First we define a triangular matrix with known eigenvalues, and a random vector $b$.

```{code-cell}
ev = 10 + arange(1, 101)
A = triu(random.rand(100, 100), 1) + diag(ev)
b = random.rand(100)
```

Next we build up the first ten Krylov matrices iteratively, using renormalization after each matrix-vector multiplication.

```{code-cell}
Km = zeros([100, 30])
Km[:, 0] = b
for m in range(29):
    v = A @ Km[:, m]
    Km[:, m + 1] = v / norm(v)
```

Now we solve least-squares problems for Krylov matrices of increasing dimension, recording the residual in each case.

```{code-cell}
from numpy.linalg import lstsq
resid = zeros(30)
resid[0] = norm(b)
for m in range(1, 30):
    z = lstsq(A @ Km[:, :m], b, rcond=None)[0]
    x = Km[:, :m] @ z
    resid[m] = norm(b - A @ x)
```

The linear system approximations show smooth linear convergence at first, but the convergence stagnates after only a few digits have been found.

```{code-cell}
semilogy(range(30), resid, "-o")
xlabel("$m$"),  ylabel("$\\| b-Ax_m \\|$")
title(("Residual for linear systems"));
```
