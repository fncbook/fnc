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
[**Demo %s**](#demo-symm-eig-rayleigh)

We will use a symmetric matrix with a known EVD and eigenvalues equal to the integers from 1 to 20.

```{code-cell}
from numpy.linalg import qr
n = 20
d = arange(n) + 1
D = diag(d)
V, _ = qr(random.randn(n, n))    # get a random orthogonal V
A = V @ D @ V.T
```

The Rayleigh quotient is a scalar-valued function of a vector.

```{code-cell}
R = lambda x: dot(x, A @ x) / dot(x, x)
```

The Rayleigh quotient evaluated at an eigenvector gives the corresponding eigenvalue.

```{code-cell}
print(R(V[:, 6]))
```

If the input to he Rayleigh quotient is within a small $\delta$ of an eigenvector, its output is within $O(\delta^2)$ of the corresponding eigenvalue. In this experiment, we observe that each additional digit of accuracy in an approximate eigenvector gives two more digits to the eigenvalue estimate coming from the Rayleigh quotient.

```{code-cell}
results = PrettyTable(["perturbation size", "R.Q. - λ"])
for delta in 1 / 10 ** arange(1, 6):
    e = random.randn(n)
    e = delta * e / norm(e)
    x = V[:, 5] + e
    quotient = R(x)
    results.add_row([delta, quotient - d[5]])

print(results)
```
