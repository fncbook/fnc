---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
Here is a previously encountered matrix that factors well.

```{code-cell}
A = [2 0 4 3 ; -4 5 -7 -10 ; 1 15 2 -4.5 ; -2 0 2 -13];
L, U = FNC.lufact(A)
L
```

If we swap the second and fourth rows of $\mathbf{A}$, the result is still nonsingular. However, the factorization now fails.

```{code-cell}
A[[2, 4], :] = A[[4, 2], :]  
L, U = FNC.lufact(A)
L
```

```{index} Julia; NaN
```

The presence of `NaN` in the result indicates that some impossible operation was required. The source of the problem is easy to locate. We can find the first outer product in the factorization just fine:

```{code-cell}
U[1, :] = A[1, :]
L[:, 1] = A[:, 1] / U[1, 1]
A -= L[:, 1] * U[1, :]'
```

The next step is `U[2, :] = A[2, :]`, which is also OK. But then we are supposed to divide by `U[2, 2]`, which is zero. The algorithm cannot continue.
