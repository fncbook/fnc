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
[**Demo %s**](#demo-structure-banded)

Here is a matrix with both lower and upper bandwidth equal to one. Such a matrix is called tridiagonal.

```{code-cell} 
A = array([ 
    [2, -1,  0,  0,  0,  0],
    [4,  2, -1,  0,  0,  0],
    [0,  3,  0, -1,  0,  0],
    [0,  0,  2,  2, -1,  0],
    [0,  0,  0,  1,  1, -1],
    [0,  0,  0,  0,  0,  2 ]
    ])
```

```{index} ! Python; diag
```

We can extract the elements on any diagonal using the `diag` command. The "main" or central diagonal is numbered zero, above and to the right of that is positive, and below and to the left is negative.

```{code-cell} 
print( diag(A) )
```

```{code-cell} 
print( diag(A, 1) )
```

```{code-cell} 
print( diag(A, -1) )
```

We can also construct matrices by specifying a diagonal with the `diag` function.

```{code-cell} 
A = A + diag([pi, 8, 6, 7], 2)
print(A)
```

```{code-cell} 
L, U = FNC.lufact(A)
print(L)
```

```{code-cell} 
print(U)
```

Observe above that the lower and upper bandwidths of $\mathbf{A}$ are preserved in the factor matrices.
