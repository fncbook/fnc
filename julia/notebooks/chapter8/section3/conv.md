---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
We set up a $5\times 5$ triangular matrix with prescribed eigenvalues on its diagonal.

```{code-cell}
λ = [1, -0.75, 0.6, -0.4, 0]
# Make a triangular matrix with eigenvalues on the diagonal.
A = triu(ones(5, 5), 1) + diagm(λ)
```

We run inverse iteration with the shift $s=0.7$ and take the final estimate as our "exact" answer to observe the convergence.

```{code-cell}
s = 0.7
β, x = FNC.inviter(A, s, 30)
eigval = β[end]
```

As expected, the eigenvalue that was found is the one closest to 0.7. The convergence is again linear.

```{code-cell}
using Plots
err = @. abs(eigval - β)
plot(0:28, err[1:end-1];
    m=:o,  xlabel=L"k", 
    yaxis=(L"|\lambda_3-\beta_k|", :log10, [1e-16, 1]),
    title="Convergence of inverse iteration")
```

The observed linear convergence rate is found from the data.

```{code-cell}
@show observed_rate = err[22] / err[21];
```

```{index} ! Julia; sortperm
```

We reorder the eigenvalues to enforce {eq}`shiftorder`.
```{tip}
The `sortperm` function returns the index permutation needed to sort the given vector, rather than the sorted vector itself.
```

```{code-cell}
λ = λ[sortperm(abs.(λ .- s))]
```

Hence the theoretical convergence rate is

```{code-cell}
@show theoretical_rate = (λ[1] - s) / (λ[2] - s);
```
