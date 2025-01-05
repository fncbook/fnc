---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
Here again is the matrix from {numref}`Demo {number} <demo-pivoting-fix>`.

```{code-cell}
A = [2 0 4 3 ; -2 0 2 -13; 1 15 2 -4.5 ; -4 5 -7 -10]
```

As the factorization proceeded, the pivots were selected from rows 4, 3, 2, and finally 1. If we were to put the rows of $\mathbf{A}$ into that order, then the algorithm would run exactly like the plain LU factorization from {numref}`section-linsys-lu`. 

```{code-cell}
B = A[[4, 3, 2, 1], :]
L, U = FNC.lufact(B);
```

We obtain the same $\mathbf{U}$ as before:

```{code-cell}
U
```

And $\mathbf{L}$ has the same rows as before, but arranged into triangular order:

```{code-cell}
L
```
