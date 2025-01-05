---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---

```{index} Julia; adjoint, Julia; \'
```

We first define a hermitian matrix. Note that the `'` operation is the adjoint and includes complex conjugation.

```{code-cell}
n = 7
A = randn(n, n) + 1im * randn(n, n)
A = (A + A') / 2
```

```{index} Julia; cond
```

We confirm that the matrix $\mathbf{A}$ is normal by checking that $\kappa(\mathbf{V}) = 1$ (to within roundoff).

```{code-cell}
λ, V = eigen(A)
@show cond(V);
```

Now we perturb $\mathbf{A}$ and measure the effect on the eigenvalues. The Bauer–Fike theorem uses absolute differences, not relative ones.
```{tip}
Since the ordering of eigenvalues can change, we look at all pairwise differences and take the minima.
```

```{code-cell}
ΔA = 1e-8 * normalize(randn(n, n) + 1im * randn(n, n))
λ̃ = eigvals(A + ΔA)
dist = minimum([abs(x - y) for x in λ̃, y in λ], dims=2)
```

As promised, the perturbations in the eigenvalues do not exceed the normwise perturbation to the original matrix.

Now we see what happens for a triangular matrix.

```{code-cell}
n = 20
x = 1:n
A = triu(x * ones(n)')
A[1:5, 1:5]
```

This matrix is not especially close to normal.

```{code-cell}
λ, V = eigen(A)
@show cond(V);
```

As a result, the eigenvalues can change by a good deal more.

```{code-cell}
ΔA = 1e-8 * normalize(randn(n, n) + 1im * randn(n, n))
λ̃ = eigvals(A + ΔA)
dist = minimum([abs(x - y) for x in λ̃, y in λ], dims=2)
BF_bound = cond(V) * norm(ΔA)
@show maximum(dist), BF_bound;
```

If we plot the eigenvalues of many perturbations, we get a cloud of points that roughly represents all the possible eigenvalues when representing this matrix with single-precision accuracy.

```{code-cell}
using Plots
plt = scatter(λ, zeros(n), aspect_ratio=1)
for _ in 1:200
    ΔA = eps(Float32) * normalize(randn(n, n) + 1im * randn(n, n))
    λ̃ = eigvals(A + ΔA)
    scatter!(real(λ̃), imag(λ̃), m=1, color=:black)
end
plt
```

The plot shows that some eigenvalues are much more affected than others. This situation is not unusual, but it is not explained by the Bauer–Fike theorem.
