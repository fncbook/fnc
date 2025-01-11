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
[**Demo %s**](#demo-pivoting-fail)

Here is a previously encountered matrix that factors well.

```{code-cell} 
A = array([
    [2, 0, 4, 3],
    [-4, 5, -7, -10],
    [1, 15, 2, -4.5],
    [-2, 0, 2, -13]
    ])
L, U = FNC.lufact(A)
print(L)
```

If we swap the second and fourth rows of $\mathbf{A}$, the result is still nonsingular. However, the factorization now fails.

```{code-cell} 
A[[1, 3], :] = A[[3, 1], :]  
L, U = FNC.lufact(A)
print(L)
```

```{index} Python; NaN
```

The presence of `NaN` in the result indicates that some impossible operation was required. The source of the problem is easy to locate. We can find the first outer product in the factorization just fine:

```{code-cell}
U[0, :] = A[0, :]
L[:, 0] = A[:, 0] / U[0, 0]
A -= outer(L[:, 0],  U[0, :])
print(A)
```

The next step is `U[1, :] = A[1, :]`, which is also OK. But then we are supposed to divide by `U[1, 1]`, which is zero. The algorithm cannot continue.
