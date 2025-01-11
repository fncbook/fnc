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
[**Demo %s**](#demo-gmres-intro)

We define a triangular matrix with known eigenvalues and a random vector $\mathbf{b}$.

```{code-cell}
ev = 10 + arange(1, 101)
A = triu(random.rand(100, 100), 1) + diag(ev)
b = random.rand(100)
```

Instead of building the Krylov matrices, we use the Arnoldi iteration to generate equivalent orthonormal vectors.

```{code-cell}
Q, H = FNC.arnoldi(A, b, 60)
print(H[:5, :5])
```

The Arnoldi bases are used to solve the least-squares problems defining the GMRES iterates.

```{code-cell}
from numpy.linalg import lstsq
resid = zeros(61)
resid[0] = norm(b)
for m in range(1, 61):
    s = hstack([norm(b), zeros(m)])
    z = lstsq(H[: m + 1, :m], s, rcond=None)[0]
    x = Q[:, :m] @ z
    resid[m] = norm(b - A @ x)
```

The approximations converge smoothly, practically all the way to machine epsilon.

```{code-cell}
semilogy(range(61), resid, "-o")
xlabel("$m$"),  ylabel("$\| b-Ax_m \|$")
title("Residual for GMRES");
```
