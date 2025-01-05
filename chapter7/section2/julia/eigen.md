---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{code-cell}
include("../../../julia/FNC_init.jl")
```
[**Demo %s**](#demo-evd-eigen)


```{index} ! Julia; eigvals
```

The `eigvals` function returns a vector of the eigenvalues of a matrix.

```{code-cell}
A = π * ones(2, 2)
```

```{code-cell}
λ = eigvals(A)
```

```{index} ! Julia; eigen
```

If you want the eigenvectors as well, use `eigen`.

```{code-cell}
λ, V = eigen(A)
```

```{code-cell}
norm(A * V[:, 2] - λ[2] * V[:, 2])
```

```{index} ! Julia; sortby
```

Both functions allow you to sort the eigenvalues by specified criteria.

```{code-cell}
A = diagm(-2.3:1.7)
@show eigvals(A, sortby=real);
@show eigvals(A, sortby=abs);
```

If the matrix is not diagonalizable, no message is given, but `V` will be singular. The robust way to detect that circumstance is via $\kappa(\mathbf{V})$.

```{index} condition number; of a matrix
```

```{code-cell}
A = [-1 1; 0 -1]
λ, V = eigen(A)
```

```{code-cell}
cond(V)
```

Even in the nondiagonalizable case, $\mathbf{A}\mathbf{V} = \mathbf{V}\mathbf{D}$ holds.

```{code-cell}
opnorm(A * V - V * diagm(λ))
```
