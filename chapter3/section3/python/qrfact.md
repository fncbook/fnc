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
[**Demo %s**](#demo-qr-qrfact)


MATLAB provides access to both the thin and full forms of the QR factorization.

```{code-cell}
A = 1.0 + floor(9 * random.rand(6,4))
A.shape
```

Here is the full form:

```{code-cell}
from numpy.linalg import qr
Q, R = qr(A, "complete")
print(f"size of Q is {Q.shape}")
print("R:")
print(R)
```

We can test that $\mathbf{Q}$ is an orthogonal matrix:

```{code-cell}
print(f"norm of (Q^T Q - I) is {norm(Q.T @ Q - eye(6)):.3e}")
```

The default for `qr`, and the one you usually want, is the thin form.

```{code-cell}
Q_hat, R_hat = qr(A)
print(f"size of Q_hat is {Q_hat.shape}")
print("R_hat:")
print(R_hat)
```

Now $\hat{\mathbf{Q}}$ cannot be an orthogonal matrix, because it is not square, but it is still ONC. Mathematically, $\hat{\mathbf{Q}}^T \hat{\mathbf{Q}}$ is a $4\times 4$ identity matrix.

```{code-cell}
print(f"norm of (Q_hat^T Q_hat - I) is {norm(Q_hat.T @ Q_hat - eye(4)):.3e}")
```
