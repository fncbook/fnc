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
[**Demo %s**](#demo-systems-backslash)

For a square matrix $A$, the command `solve(A, B)` from `numpy.linalg` is mathematically equivalent to $\mathbf{A}^{-1} \mathbf{b}$. 

```{code-cell} 
A = array([[1, 0, -1], [2, 2, 1], [-1, -3, 0]])
b = array([1, 2, 3])
```

```{code-cell}
x = linalg.solve(A, b)
print(x)
```

```{index} residual
```

One way to check the answer is to compute a quantity known as the **residual**. It is (ideally) close to machine precision(relative to the elements in the data). 

```{code-cell} 
residual = b - A @ x
print(residual)
```

If the matrix $\mathbf{A}$ is singular, you may get an error.

```{code-cell} 
:tags: [raises-exception]
A = array([[0, 1], [0, 0]])
b = array([1, -1])
linalg.solve(A, b)    # error, singular matrix
```

A linear system with a singular matrix might have no solution or infinitely many solutions, but in either case, a numerical solution becomes trickier. Detecting singularity is a lot like checking whether two floating-point numbers are *exactly* equal: because of roundoff, it could be missed. We're headed toward a more robust way to fully describe this situation.
