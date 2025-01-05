---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
A system of nonlinear equations is defined by its residual and Jacobian.
```{tip}
Be careful when coding a Jacobian all in one statement. Spaces separate columns, so `x[3]-1` is not the same as `x[3] - 1`.
```

```{code-cell}
function func(x)
    [exp(x[2] - x[1]) - 2,
        x[1] * x[2] + x[3],
        x[2] * x[3] + x[1]^2 - x[2]
    ]
end;

function jac(x)
    [
        -exp(x[2] - x[1]) exp(x[2] - x[1]) 0
        x[2] x[1] 1
        2*x[1] x[3]-1 x[2]
    ]
end;
```

We will use a `BigFloat` starting value, and commensurately small stopping tolerances, in order to get a sequence long enough to measure convergence.

```{code-cell}
x₁ = BigFloat.([0, 0, 0])
ϵ = eps(BigFloat)
x = FNC.newtonsys(func, jac, x₁, xtol=ϵ, ftol=ϵ);
```

Let's compute the residual of the last result in order to check the quality.

```{code-cell}
r = x[end]
@show residual = norm(func(r));
```

We take the sequence of norms of errors, applying the log so that we can look at the exponents.

```{code-cell}
logerr = [Float64(log(norm(r - x[k]))) for k in 1:length(x)-1]
ratios = [NaN; [logerr[i+1] / logerr[i] for i in 1:length(logerr)-1]]
@pt :header=["iteration", "log error", "ratio"] [eachindex(logerr) logerr ratios]
```

The ratio is neatly converging toward 2, which is expected for quadratic convergence.
