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
[**Demo %s**](#demo-evd-bauerfike)


We first define a hermitian matrix. Note that we add the *conjugate* transpose of a matrix to itself.

```{code-cell}
n = 7
A = random.randn(n, n) + 1j * random.randn(n, n)
A = (A + conj(A.T)) / 2
```

```{index} Python; cond
```

We confirm that the matrix $\mathbf{A}$ is normal by checking that $\kappa(\mathbf{V}) = 1$ (to within roundoff).

```{code-cell}
from numpy.linalg import eig, cond
d, V = eig(A)
print(f"eigenvector matrix has condition number {cond(V):.5f}")
```

Now we perturb $\mathbf{A}$ and measure the effect on the eigenvalues. Note that the Bauer–Fike theorem uses absolute differences, not relative ones. Since the ordering of eigenvalues can change, we look at all pairwise differences and take the minima.

```{code-cell}
E = random.randn(n, n) + 1j * random.randn(n, n)
E = 1e-8 * E / norm(E, 2)
dd, _ = eig(A + E)
dist = array([min([abs(x - y) for x in dd]) for y in d])
print(dist)
```
As promised, the perturbations in the eigenvalues do not exceed the normwise perturbation to the original matrix.

Now we see what happens for a triangular matrix.

```{code-cell}
n = 20
x = arange(n) + 1
A = triu(outer(x, ones(n)))
print(A[:5, :5])
```

This matrix is not at all close to normal.

```{code-cell}
d, V = eig(A)
print(f"eigenvector matrix has condition number {cond(V):.2e}")
```

As a result, the eigenvalues can change by a good deal more.

```{code-cell}
E = random.randn(n, n) + 1j * random.randn(n, n)
E = 1e-8 * E / norm(E, 2)
dd, _ = eig(A + E)
dist = array([min([abs(x - y) for x in dd]) for y in d])
print(f"Maximum eigenvalue change is {max(dist):.2e}")
print(f"The Bauer-Fike upper bound is {cond(V) * norm(E, 2):.2e}")
```

If we plot the eigenvalues of many perturbations, we get a cloud of points that roughly represents all the possible eigenvalues when representing this matrix with single-precision accuracy.

```{code-cell}
clf
scatter(d, zeros(n), 18)
axis("equal") 
for _ in range(100):
    E = random.randn(n, n) + 1j * random.randn(n, n)
    E = finfo(np.float32).eps * E / norm(E, 2)
    dd, _ = eig(A + E)
    scatter(real(dd), imag(dd), 2, 'k')
```

The plot shows that some eigenvalues are much more affected than others. This situation is not unusual, but it is not explained by the Bauer–Fike theorem.
