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
[**Demo %s**](#demo-structure-linalg)


The following generates a random sparse matrix with prescribed eigenvalues.

```{code-cell}
n = 4000
density = 4e-4
ev = 1 / arange(1, n + 1)
A = FNC.sprandsym(n, density, eigvals=ev)
print(f"density is {A.nnz / prod(A.shape):.3%}")
```

```{index} ! Python; eigs
```

The `eigs` function finds a small number eigenvalues meeting some criterion. First, we ask for the 5 of largest (complex) magnitude using `which="LM"`.

```{code-cell}
from scipy.sparse.linalg import eigs
ev, V = eigs(A, k=5, which="LM")    # largest magnitude
print(1 / ev)
```

Now we find the 4 closest to the value 1 in the complex plane, via `sigma=1`.

```{code-cell}
from scipy.sparse.linalg import eigs
ev, V = eigs(A, k=4, sigma=0.03)    # closest to sigma
print(ev)
```

The time needed to solve a sparse linear system is not easy to predict unless you have some more information about the matrix. But it will typically be orders of magnitude faster than the dense version of the same problem.

```{code-cell}
from scipy.sparse.linalg import spsolve
x = 1 / arange(1, n + 1)
b = A @ x
start = timer()
xx = spsolve(A, b)
print(f"sparse time: {timer() - start:.3g} sec")
print(f"residual: {norm(b - A @ xx, 2):.1e}")
```

```{code-cell}
from numpy.linalg import solve
F = A.todense()
start = timer()
xx = solve(F, b)
print(f"dense time: {timer() - start:.3g} sec")
print(f"residual: {norm(b - A @ xx, 2):.1e}")
```
