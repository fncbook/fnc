---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
The third output of `plufact` is the permutation vector we need to apply to $\mathbf{A}$.

```{code-cell}
A = rand(1:20, 4, 4)
L, U, p = FNC.plufact(A)
A[p,:] - L * U   # should be ≈ 0
```

Given a vector $\mathbf{b}$, we solve $\mathbf{A}\mathbf{x}=\mathbf{b}$ by first permuting the entries of $\mathbf{b}$ and then proceeding as before.

```{code-cell}
b = rand(4)
z = FNC.forwardsub(L,b[p])
x = FNC.backsub(U,z)
```

A residual check is successful:

```{code-cell}
b - A*x
```
