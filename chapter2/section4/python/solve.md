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
[**Demo %s**](#demo-lu-solve)

Here are the data for a linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$. 

```{code-cell}
A = array([
    [2, 0, 4, 3], 
    [-4, 5, -7, -10], 
    [1, 15, 2, -4.5],
    [-2, 0, 2, -13]
    ])
b = array([4, 9, 9, 4])
```

We apply {numref}`Function {number} <function-lufact>` and then do two triangular solves.

```{code-cell}
L, U = FNC.lufact(A)
z = FNC.forwardsub(L, b)
x = FNC.backsub(U, z)
```

A check on the residual assures us that we found the solution.

```{code-cell}
b - A @ x
```

