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
[**Demo %s**](#demo-pivoting-stable)

We construct a linear system for this matrix with $\epsilon=10^{-12}$ and exact solution $[1,1]$:

```{code-cell}
ep = 1e-12
A = array([[-ep, 1], [1, -1]])
b = A @ array([1, 1])
```

We can factor the matrix without pivoting and solve for $\mathbf{x}$.

```{code-cell}
L, U = FNC.lufact(A)
print(FNC.backsub( U, FNC.forwardsub(L, b) ))
```

Note that we have obtained only about 5 accurate digits for $x_1$. We could make the result even more inaccurate by making $\epsilon$ even smaller:

```{code-cell}
ep = 1e-20;
A = array([[-ep, 1], [1, -1]])
b = A @ array([1, 1])
L, U = FNC.lufact(A)
print(FNC.backsub( U, FNC.forwardsub(L, b) ))
```

This effect is not due to ill conditioning of the problem—a solution with PLU factorization works perfectly:

```{code-cell}
print(linalg.solve(A, b))
```
