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
[**Demo %s**](#demo-pivoting-builtin)

In `linalg.solve`, the matrix `A` is PLU-factored, followed by two triangular solves. If we want to do those steps seamlessly, we can use the `lu_factor` and `lu_solve` from `scipy.linalg`.

```{code-cell}
from scipy.linalg import lu_factor, lu_solve
A = random.randn(500, 500) 
b = ones(500)  
LU, perm = lu_factor(A)
x = lu_solve((LU, perm), b)
```

Why would we ever bother with this? In {numref}`section-linsys-efficiency` we showed that the factorization is by far the most costly part of the solution process. A factorization object allows us to do that costly step only once per matrix, but solve with multiple right-hand sides.

```{code-cell}
start = timer()
for k in range(50): linalg.solve(A, random.rand(500))
print(f"elapsed time for 50 full solves: {timer() - start}")

start = timer()
LU, perm = lu_factor(A)
for k in range(50): lu_solve((LU, perm), random.rand(500))
print(f"elapsed time for 50 shortcut solves: {timer() - start}")
```
