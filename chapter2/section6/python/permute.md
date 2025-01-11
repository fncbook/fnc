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
[**Demo %s**](#demo-pivoting-permute)

Here again is the matrix from {numref}`Demo {number} <demo-pivoting-fix>`.

```{code-cell}
A = array([
    [2, 0, 4, 3],
    [-2, 0, 2, -13],
    [1, 15, 2, -4.5],
    [-4, 5, -7, -10]
    ])
```

As the factorization proceeded, the pivots were selected from rows 4, 3, 2, and finally 1 (with NumPy indices being one less). If we were to put the rows of $\mathbf{A}$ into that order, then the algorithm would run exactly like the plain LU factorization from {numref}`section-linsys-lu`. 

```{code-cell}
B = A[[3, 2, 1, 0], :]
L, U = FNC.lufact(B);
```

We obtain the same $\mathbf{U}$ as before:

```{code-cell}
print(U)
```

And $\mathbf{L}$ has the same rows as before, but arranged into triangular order:

```{code-cell}
print(L)
```
