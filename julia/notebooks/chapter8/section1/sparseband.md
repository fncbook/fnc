---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---

```{index} ! Julia; spdiagm
```

The `spdiagm` function creates a sparse matrix given its diagonal elements. The main or central diagonal is numbered zero, above and to the right of that is positive, and below and to the left is negative.

```{code-cell}
n = 50;
A = spdiagm(-3 => fill(n, n - 3),
    0 => ones(n),
    1 => -(1:n-1),
    5 => fill(0.1, n - 5))
Matrix(A[1:7, 1:7])
```

```{index} ! Julia; sparse
```


:::{grid-item}


Without pivoting, the LU factors have the same lower and upper bandwidth as the original matrix.


:::
```{tip}
The `sparse` function converts any matrix to sparse form. But it's usually better to construct a sparse matrix directly, as the standard form might not fit in memory.
```

```{code-cell}
L, U = FNC.lufact(A)
plot(layout=2)
spy!(sparse(L), m=2, subplot=1, title=L"\mathbf{L}", color=:blues)
spy!(sparse(U), m=2, subplot=2, title=L"\mathbf{U}", color=:blues)
```

However, if we introduce row pivoting, bandedness may be expanded or destroyed.

```{code-cell}
fact = lu(A)
plot(layout=2)
spy!(sparse(fact.L), m=2, subplot=1, title=L"\mathbf{L}", color=:blues)
spy!(sparse(fact.U), m=2, subplot=2, title=L"\mathbf{U}", color=:blues)
```
