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
[**Demo %s**](#demo-structure-fill)


Here is the adjacency matrix of a graph representing a small-world network, featuring connections to neighbors and a small number of distant contacts.

```{code-cell}
import networkx as nx
wsg = nx.watts_strogatz_graph(200, 4, 0.02)
```

Because each node connects to relatively few others, the adjacency matrix is quite sparse.

```{code-cell}
A = nx.adjacency_matrix(wsg)
spy(A)
title("Adjacency matrix $A$");
```

By {numref}`Theorem {number} <theorem-insight-adjmat>`, the entries of $\mathbf{A}^k$ give the number of walks of length $k$ between pairs of nodes, as with "*k* degrees of separation" within a social network. As $k$ grows, the density of $\mathbf{A}^k$ also grows.
```{tip}
:class: dropdown
While `A**6` is valid syntax here, it means elementwise power, not matrix power. 
```

```{index} ! Python; matrix_power
```

```{code-cell}
from scipy.sparse.linalg import matrix_power
spy(matrix_power(A, 6))
title(("$A^6$"));
```
