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
[**Demo %s**](#demo-evd-eigen)


```{index} ! Python; eig
```

The `eig` function from `scipy.linalg` will return a vector of eigenvalues and a matrix of associated eigenvectors.

```{code-cell}
from numpy.linalg import eig
A = pi * ones([2, 2])
d, V = eig(A)
print("eigenvalues:", d)
```

We can check the fact that this is an EVD (although in practice we never invert a matrix).

```{code-cell}
from numpy.linalg import inv
D = diag(d)
print(f"should be near zero: {norm(A - V @ D @ inv(V), 2):.2e}")
```

If the matrix is not diagonalizable, no message is given, but `V` will be singular. The robust way to detect that circumstance is via $\kappa(\mathbf{V})$.

```{index} condition number; of a matrix
```

```{index} Python; cond
```

```{code-cell}
from numpy.linalg import cond
A = array([[1, 1], [0, 1]])
d, V = eig(A)
print(f"cond(V) is {cond(V):.2e}")
```

But even in the nondiagonalizable case, $\mathbf{A}\mathbf{V} = \mathbf{V}\mathbf{D}$ holds up to roundoff error.

```{code-cell}
print(f"should be near zero: {norm(A @ V - V @ diag(d), 2):.2e}")
```
```
