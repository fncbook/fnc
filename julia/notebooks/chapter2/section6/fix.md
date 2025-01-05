---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
Here is the trouble-making matrix from {numref}`Demo {number} <demo-pivoting-fail>`.

```{code-cell}
A₁ = [2 0 4 3 ; -2 0 2 -13; 1 15 2 -4.5 ; -4 5 -7 -10]
```

We now find the largest candidate pivot in the first column. We don't care about sign, so we take absolute values before finding the max.
```{tip}
The `argmax` function returns the location of the largest element of a vector or matrix.
```


```{code-cell}
i = argmax( abs.(A₁[:, 1]) ) 
```

This is the row of the matrix that we extract to put into $\mathbf{U}$. That guarantees that the division used to find $\boldsymbol{\ell}_1$ will be valid.

```{code-cell}
L, U = zeros(4,4),zeros(4,4)
U[1, :] = A₁[i, :]
L[:, 1] = A₁[:, 1] / U[1, 1]
A₂ = A₁ - L[:, 1] * U[1, :]'
```

Observe that $\mathbf{A}_2$ has a new zero row and zero column, but the zero row is the fourth rather than the first. However, we forge on by using the largest possible pivot in column 2 for the next outer product.

```{code-cell}
@show i = argmax( abs.(A₂[:, 2]) ) 
U[2, :] = A₂[i, :]
L[:, 2] = A₂[:, 2] / U[2, 2]
A₃ = A₂ - L[:, 2] * U[2, :]'
```

Now we have zeroed out the third row as well as the second column. We can finish out the procedure.

```{code-cell}
@show i = argmax( abs.(A₃[:, 3]) ) 
U[3, :] = A₃[i, :]
L[:, 3] = A₃[:, 3] / U[3, 3]
A₄ = A₃ - L[:, 3] * U[3, :]'
```

```{code-cell}
@show i = argmax( abs.(A₄[:, 4]) ) 
U[4, :] = A₄[i, :]
L[:, 4] = A₄[:, 4] / U[4, 4];
```

We do have a factorization of the original matrix:

```{code-cell}
A₁ - L * U
```

And $\mathbf{U}$ has the required structure:

```{code-cell}
U
```

However, the triangularity of $\mathbf{L}$ has been broken.

```{code-cell}
L
```
