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
[**Demo %s**](#demo-structure-cholesky)

A randomly chosen matrix is extremely unlikely to be symmetric. However, there is a simple way to symmetrize one.

```{code-cell} 
A = 1.0 + floor(9 * random.rand(4, 4))
B = A + A.T
print(B)
```

Similarly, a random symmetric matrix is unlikely to be positive definite. The Cholesky algorithm always detects a non-PD matrix by quitting with an error.

```{index} ! Python; cholesky
```

```{code-cell} 
:tags: [raises-exception]
from numpy.linalg import cholesky
cholesky(B)    # raises an exception, not positive definite
```

It's not hard to manufacture an SPD matrix to try out the Cholesky factorization:

```{code-cell} 
B = A.T @ A
R = cholesky(B)
print(R)
```

```{code-cell} 
print(norm(R @ R.T - B) / norm(B))
```
