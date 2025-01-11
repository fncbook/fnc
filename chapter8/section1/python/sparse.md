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
[**Demo %s**](#demo-structure-sparse)

```{tip}
Functions to work with sparse matrices are found in the `scipy.sparse` module.
```

Here we load the adjacency matrix of a graph with 2790 nodes. Each node is a web page referring to Roswell, NM, and the edges represent links between web pages. (Credit goes to Panayiotis Tsaparas and the University of Toronto for making this data public.)

```{code-cell}
import scipy.sparse as sp
from scipy.io import loadmat

vars = loadmat("roswelladj.mat")    # get from the book's website
A = sp.csr_matrix(vars["A"])
```

We may define the density of $\mathbf{A}$ as the number of nonzeros divided by the total number of entries.

```{index} ! Python; nnz
```

```{code-cell}
m, n = A.shape
print(f"density is {A.nnz / (m * n):.3%}")
```

We can compare the storage space needed for the sparse $\mathbf{A}$ with the space needed for its dense / full counterpart.


```{code-cell}
F = A.todense()
print(f"{A.data.nbytes/1e6:.3f} MB for sparse form, {F.nbytes/1e6:.3f} MB for dense form")
```

Matrix-vector products are also much faster using the sparse form because operations with structural zeros are skipped.

```{code-cell}
from timeit import default_timer as timer
x = random.randn(n)
start = timer()
for i in range(1000):
    A @ x
print(f"sparse time: {timer() - start:.4g} sec")
```

```{code-cell}
start = timer()
for i in range(1000):
    F @ x
print(f"dense time: {timer() - start:.4g} sec")
```
