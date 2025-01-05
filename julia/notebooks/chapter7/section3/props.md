---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
We verify some of the fundamental SVD properties using standard Julia functions from `LinearAlgebra`.

```{code-cell}
A = [i^j for i = 1:5, j = 0:3]
```

```{index} ! Julia; svdvals
```


To get only the singular values, use `svdvals`.

```{code-cell}
σ = svdvals(A)
```

Here is verification of the connections between the singular values, norm, and condition number.

```{code-cell}
@show opnorm(A, 2);
@show σ[1];
```

```{code-cell}
@show cond(A, 2);
@show σ[1] / σ[end];
```

```{index} ! Julia; svd
```

To get singular vectors as well, use `svd`. The thin form of the factorization is the default.

```{code-cell}
U, σ, V = svd(A);
@show size(U);
@show size(V);
```

We verify the orthogonality of the singular vectors as follows:

```{code-cell}
@show opnorm(U' * U - I);
@show opnorm(V' * V - I);
```
