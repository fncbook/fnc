---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
A randomly chosen matrix is extremely unlikely to be symmetric. However, there is a simple way to symmetrize one.

```{code-cell}
A = rand(1.0:9.0, 4, 4)
B = A + A'
```

```{index} ! Julia; cholesky
```

Similarly, a random symmetric matrix is unlikely to be positive definite. The Cholesky algorithm always detects a non-PD matrix by quitting with an error.
```{tip}
The `cholesky` function computes a Cholesky factorization if possible, or throws an error for a non-positive-definite matrix. However, it does *not* check for symmetry.
```

```{code-cell}
:tags: [raises-exception]
cholesky(B)    # throws an error
```

It's not hard to manufacture an SPD matrix to try out the Cholesky factorization.

```{code-cell}
B = A' * A
cf = cholesky(B)
```

What's returned is a factorization object. Another step is required to extract the factor as a matrix:

```{code-cell}
R = cf.U
```

Here we validate the factorization:

```{code-cell}
opnorm(R' * R - B) / opnorm(B)
```
