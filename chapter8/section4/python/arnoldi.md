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
[**Demo %s**](#demo-subspace-arnoldi)

We illustrate a few steps of the Arnoldi iteration for a small matrix.

```{code-cell}
A = random.choice(range(10), (6, 6))
print(A)
```

The seed vector we choose here determines the first member of the orthonormal basis.

```{code-cell}
u = random.randn(6)
Q = zeros([6, 3])
Q[:, 0] = u / norm(u)
```

Multiplication by $\mathbf{A}$ gives us a new vector in $\mathcal{K}_2$.

```{code-cell}
Aq = A @ Q[:, 0]
```

We subtract off its projection in the previous direction. The remainder is rescaled to give us the next orthonormal column.

```{code-cell}
v = Aq - dot(Q[:, 0], Aq) * Q[:, 0]
Q[:, 1] = v / norm(v)
```

On the next pass, we have to subtract off the projections in two previous directions.

```{code-cell}
Aq = A @ Q[:, 1]
v = Aq - dot(Q[:, 0], Aq) * Q[:, 0] - dot(Q[:, 1], Aq) * Q[:, 1]
Q[:, 2] = v / norm(v)
```

At every step, $\mathbf{Q}_m$ is an ONC matrix.

```{code-cell}
print(f"should be near zero: {norm(Q.T @ Q - eye(3)):.2e}")
```

And $\mathbf{Q}_m$ spans the same space as the three-dimensional Krylov matrix.

```{code-cell}
from numpy.linalg import matrix_rank
K = stack([u, A @ u, A @ A @ u], axis=-1)
Q_and_K = hstack([Q, K])
print(matrix_rank(Q_and_K))
```
