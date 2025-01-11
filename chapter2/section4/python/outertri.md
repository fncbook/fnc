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
[**Demo %s**](#demo-lu-outertri)

```{index} Python; tril, Python; triu
```
We explore the outer product formula for two random triangular matrices.

```{code-cell}
from numpy.random import randint
L = tril(randint(1, 10, size=(3, 3)))
print(L)
```

```{code-cell}
U = triu(randint(1, 10, size=(3, 3)))
print(U)
```

Here are the three outer products appearing in the sum in {eq}`matrixouter`:

```{code-cell}
print(outer(L[:, 0], U[0, :]))
```

```{code-cell}
print(outer(L[:, 1], U[1, :]))
```

```{code-cell}
print(outer(L[:, 2], U[2, :]))
```

Simply because of the triangular zero structures, only the first outer product contributes to the first row and first column of the entire product. 
