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
[**Demo %s**](#demo-precond-gmres)

Here is a random nonsymmetric matrix.

```{code-cell}
import scipy.sparse as sp
n = 8000
A = 2.8 * sp.eye(n) + sp.rand(n, n, 0.002)
```

Without a preconditioner, GMRES can solve a system with this matrix.

```{code-cell}
from scipy.sparse.linalg import gmres

b = random.rand(n)
hist = lambda rvec: resid.append(norm(rvec))
resid = [1.]

start = timer()
x, flag = gmres(A, b, maxiter=300, rtol=1e-10, restart=50, callback=hist)
print(f"time for plain GMRES: {timer() - start:.3f} sec")
resid_plain = resid.copy()
```

```{index} ! Python; spilu
```

The following version of incomplete LU factorization drops all sufficiently small values (i.e., replaces them with zeros). This keeps the number of nonzeros in the factors under control.

```{code-cell}
from scipy.sparse.linalg import spilu
iLU = spilu(A, drop_tol=0.2)
print(f"Factors have {iLU.nnz} nonzeros, while A has {A.nnz}")
```

The result is not a true factorization of the original matrix. However, it's close enough for an approximate inverse in a preconditioner. 

```{code-cell}
from scipy.sparse.linalg import LinearOperator
prec = LinearOperator((n, n), matvec=lambda y: iLU.solve(y))

resid = [1.];  start = timer()
x, flag = gmres(A, b, M=prec, maxiter=300, rtol=1e-10, restart=50, callback=hist)
print(f"time for preconditioned GMRES: {timer() - start:.3f} sec")
resid_prec = resid
```

```{code-cell}
:tags: [hide-input]
semilogy(resid_plain, label="no prec.")
semilogy(resid_prec, label="iLU prec.")
xlabel("iteration number"),  ylabel("residual norm")
legend()
title("GMRES convergence compared");
``` 
