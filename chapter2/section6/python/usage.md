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
[**Demo %s**](#demo-pivoting-usage)

The third output of `plufact` is the permutation vector we need to apply to $\mathbf{A}$.

```{code-cell}
A = random.randn(4, 4)
L, U, p = FNC.plufact(A)
A[p, :] - L @ U   # should be ≈ 0
```

Given a vector $\mathbf{b}$, we solve $\mathbf{A}\mathbf{x}=\mathbf{b}$ by first permuting the entries of $\mathbf{b}$ and then proceeding as before.

```{code-cell}
b = random.randn(4)
z = FNC.forwardsub(L, b[p])
x = FNC.backsub(U, z)
```

A residual check is successful:

```{code-cell}
b - A @ x
```
