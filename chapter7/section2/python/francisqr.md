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
[**Demo %s**](#demo-evd-francisqr)

Let's start with a known set of eigenvalues and an orthogonal eigenvector basis.

```{code-cell}
from numpy.linalg import qr
D = diag([-6, -1, 2, 4, 5])
V, R = qr(random.randn(5, 5))
A = V @ D @ V.T    # note that V.T = inv(V) here
```

```{code-cell}
print(sort(eig(A)[0]))
```

Now we will take the QR factorization and just reverse the factors.

```{code-cell}
Q, R = qr(A)
A = R @ Q;
```

It turns out that this is a similarity transformation, so the eigenvalues are unchanged.

```{code-cell}
print(sort(eig(A)[0]))
```

What's remarkable, and not elementary, is that if we repeat this transformation many times, the resulting matrix converges to $\mathbf{D}$.

```{code-cell}
for k in range(40):
    Q, R = qr(A)
    A = R @ Q
set_printoptions(precision=4)
print(A)
```
