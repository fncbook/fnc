---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
We will experiment with the power iteration on a 5×5 matrix with prescribed eigenvalues and dominant eigenvalue at 1.

```{code-cell}
λ = [1, -0.75, 0.6, -0.4, 0]
# Make a triangular matrix with eigenvalues on the diagonal.
A = triu(ones(5, 5), 1) + diagm(λ)
```

We run the power iteration 60 times. The best estimate of the dominant eigenvalue is the last entry of the first output.

```{code-cell}
β, x = FNC.poweriter(A, 60)
eigval = β[end]
```

We check for linear convergence using a log-linear plot of the error.

```{code-cell}
using Plots
err = @. 1 - β
plot(0:59, abs.(err); m=:o, 
    xlabel=L"k",  
    yaxis=(L"|\lambda_1-\beta_k|", :log10, [1e-10, 1]),
    title="Convergence of power iteration")
```

The asymptotic trend seems to be a straight line, consistent with linear convergence. To estimate the convergence rate, we look at the ratio of two consecutive errors in the linear part of the convergence curve. The ratio of the first two eigenvalues should match the observed rate.

```{code-cell}
@show theory = λ[2] / λ[1];
@show observed = err[40] / err[39];
```

Note that the error is supposed to change sign on each iteration. The effect of these alternating signs is that estimates oscillate around the exact value.

```{code-cell}
β[26:30]
```

In practical situations, we don't know the exact eigenvalue that the algorithm is supposed to find. In that case we would base errors on the final $\beta$ that was found, as in the following plot.

```{code-cell}
err = @. β[end] - β[1:end-1]
plot(0:58, abs.(err), m=:o, 
    xlabel=L"k", 
    yaxis=(L"|\beta_{60}-\beta_k|", :log10, [1e-10, 1]),
    title="Convergence of power iteration")
```

The results are very similar until the last few iterations, when the limited accuracy of the reference value begins to show. That is, while it is a good estimate of $\lambda_1$, it is less good as an estimate of the error in nearby estimates.

