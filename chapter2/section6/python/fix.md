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
[**Demo %s**](#demo-pivoting-fix)

Here is the trouble-making matrix from {numref}`Demo {number} <demo-pivoting-fail>`.

```{code-cell}
A_1 = array([
    [2, 0, 4, 3],
    [-2, 0, 2, -13],
    [1, 15, 2, -4.5],
    [-4, 5, -7, -10]
    ])
```

We now find the largest candidate pivot in the first column. We don't care about sign, so we take absolute values before finding the max.
```{tip}
:class: dropdown
The `argmax` function returns the location of the largest element of a vector or matrix.
```


```{code-cell}
i = argmax( abs(A_1[:, 0]) )
print(i)
```

This is the row of the matrix that we extract to put into $\mathbf{U}$. That guarantees that the division used to find $\boldsymbol{\ell}_1$ will be valid.

```{code-cell}
L, U = eye(4), zeros((4, 4))
U[0, :] = A_1[i, :]
L[:, 0] = A_1[:, 0] / U[0, 0]
A_2 = A_1 - outer(L[:, 0], U[0, :])
print(A_2)
```

Observe that $\mathbf{A}_2$ has a new zero row and zero column, but the zero row is the fourth rather than the first. However, we forge on by using the largest possible pivot in column 2 for the next outer product.

```{code-cell}
i = argmax( abs(A_2[:, 1]) ) 
print(f"new pivot row is {i}")
U[1, :] = A_2[i, :]
L[:, 1] = A_2[:, 1] / U[1, 1]
A_3 = A_2 - outer(L[:, 1], U[1, :])
print(A_3)
```
Now we have zeroed out the third row as well as the second column. We can finish out the procedure.

```{code-cell}
i = argmax( abs(A_3[:, 2]) ) 
print(f"new pivot row is {i}")
U[2, :] = A_3[i, :]
L[:, 2] = A_3[:, 2] / U[2, 2]
A_4 = A_3 - outer(L[:, 2], U[2, :])
print(A_4)
```

```{code-cell}
i = argmax( abs(A_4[:, 3]) ) 
print(f"new pivot row is {i}")
U[3, :] = A_4[i, :]
L[:, 3] = A_4[:, 3] / U[3, 3];
```

We do have a factorization of the original matrix:

```{code-cell}
A_1 - L @ U
```

And $\mathbf{U}$ has the required structure:

```{code-cell}
print(U)
```

However, the triangularity of $\mathbf{L}$ has been broken.

```{code-cell}
print(L)
```
