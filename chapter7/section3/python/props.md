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
[**Demo %s**](#demo-svd-props)

We verify some of the fundamental SVD properties using standard Julia functions from `LinearAlgebra`.

```{code-cell}
A = array([[(i + 1.0) ** j for j in range(4)] for i in range(5)])
set_printoptions(precision=4)
print(A)
```

```{index} ! Python; svd
```

The factorization is obtained using `svd` from `numpy.linalg`.

```{code-cell}
from numpy.linalg import svd
U, sigma, Vh = svd(A)
print("singular values:")
print(sigma)
```

By default, the full factorization type is returned. This can be a memory hog if one of the dimensions of $\mathbf{A}$ is very large.

```{code-cell}
print("size of U:", U.shape)
print("size of V:", Vh.T.shape)
```

Both $\mathbf{U}$ and $\mathbf{V}$ are orthogonal (in the complex case, unitary). Note that it's $\mathbf{V}^*$ that is returned, not $\mathbf{V}$.

```{code-cell}
print(f"should be near zero: {norm(U.T @ U - eye(5), 2):.2e}")
print(f"should be near zero: {norm(Vh @ Vh.T - eye(4), 2):.2e}")
```

Next we test that we have the factorization promised by the SVD, using `diagsvd` to construct a rectangular diagonal matrix.

```{code-cell}
from scipy.linalg import diagsvd
S = diagsvd(sigma, 5, 4)
print(f"should be near zero: {norm(A - U @ S @ Vh, 2):.2e}")
```

Here is verification of the connections between the singular values, norm, and condition number.

```{code-cell}
from numpy.linalg import cond
print("largest singular value:", sigma[0])
print("2-norm of the matrix:  ", norm(A, 2))
print("singular value ratio:", sigma[0] / sigma[-1])
print("2-norm condition no.:", cond(A, 2))
```

For matrices that are much taller than they are wide, the thin SVD form is more memory-efficient, because $\mathbf{U}$ takes the same shape.

```{code-cell}
A = random.randn(1000, 10)
U, sigma, Vh = svd(A, full_matrices=False)
print("size of U:", U.shape)
print("size of V:", Vh.shape)
```
