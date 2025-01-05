---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{index} ! Julia; fill, Julia; diagm, ! Julia; diag
```

Here is a small tridiagonal matrix. Note that there is one fewer element on the super- and subdiagonals than on the main diagonal.
```{tip}
Use `fill` to create an array of a given size, with each element equal to a provided value.
```

```{code-cell}
A = diagm( -1 => [4, 3, 2, 1, 0], 
    0 => [2, 2, 0, 2, 1, 2], 
    1 => fill(-1, 5) )
```

```{index} ! Julia; diag
```

We can extract the elements on any diagonal using the `diag` function. The main or central diagonal is numbered zero, above and to the right of that is positive, and below and to the left is negative.
```{tip}
The `diag` function extracts the elements from a specified diagonal of a matrix.
```

```{code-cell}
@show diag_main = diag(A);
@show diag_minusone = diag(A, -1);
```
The lower and upper bandwidths of $\mathbf{A}$ are repeated in the factors from the unpivoted LU factorization. 

```{code-cell}
L, U = FNC.lufact(A)
L
```

```{code-cell}
U
```
