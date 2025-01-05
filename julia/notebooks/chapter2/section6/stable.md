---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
We construct a linear system for this matrix with $\epsilon=10^{-12}$ and exact solution $[1,1]$:

```{code-cell}
ϵ = 1e-12
A = [-ϵ 1; 1 -1]
b = A * [1, 1]
```

We can factor the matrix without pivoting and solve for $\mathbf{x}$.

```{code-cell}
L, U = FNC.lufact(A)
x = FNC.backsub( U, FNC.forwardsub(L, b) )
```

Note that we have obtained only about 5 accurate digits for $x_1$. We could make the result even more inaccurate by making $\epsilon$ even smaller:

```{code-cell}
ϵ = 1e-20; A = [-ϵ 1; 1 -1]
b = A * [1, 1]
L, U = FNC.lufact(A)
x = FNC.backsub( U, FNC.forwardsub(L, b) )
```

This effect is not due to ill conditioning of the problem—a solution with PLU factorization works perfectly:

```{code-cell}
A \ b
```
