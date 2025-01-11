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
[**Demo %s**](#demo-structure-symm)

We begin with a symmetric $\mathbf{A}$. 

```{code-cell} 
A_1 = array([
    [2,     4,     4,     2],
    [4,     5,     8,    -5],
    [4,     8,     6,     2],
    [2,    -5,     2,   -26]
    ])
```

We won't use pivoting, so the pivot element is at position (1,1). This will become the first element on the diagonal of $\mathbf{D}$. Then we divide by that pivot to get the first column of $\mathbf{L}$.

```{code-cell}
L = eye(4)
d = zeros(4)
d[0] = A_1[0, 0]
L[:, 0] = A_1[:, 0] / d[0]
A_2 = A_1 - d[0] * outer(L[:, 0], L[:, 0])
print(A_2)
```
We are now set up the same way for the submatrix in rows and columns 2–4.

```{code-cell}
d[1] = A_2[1, 1]
L[:, 1] = A_2[:, 1] / d[1]
A_3 = A_2 - d[1] * outer(L[:, 1], L[:, 1])
print(A_3)
```

We continue working our way down the diagonal.

```{code-cell}
d[2] = A_3[2, 2]
L[:, 2] = A_3[:, 2] / d[2]
A_4 = A_3 - d[2] * outer(L[:, 2], L[:, 2])
print(A_4)
```

We have arrived at the desired factorization.

```{code-cell}
d[3] = A_4[3, 3]
print("diagonal of D:")
print(d)
print("L:")
print(L)
```

This should be comparable to machine roundoff:

```{code-cell}
print(norm(A_1 - (L @ diag(d) @ L.T), 2) / norm(A_1))
```
