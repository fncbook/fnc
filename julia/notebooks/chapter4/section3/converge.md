---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
We again look at finding a solution of $x e^x=2$ near $x=1$. To apply Newton's method, we need to calculate values of both the residual function $f$ and its derivative.

```{code-cell}
f(x) = x * exp(x) - 2;
df_dx(x) = exp(x) * (x + 1);
```

We don't know the exact root, so we use `nlsolve` to determine a proxy for it.

```{code-cell}
using NLsolve
r = nlsolve(x -> f(x[1]), [1.0]).zero
```

We use $x_1=1$ as a starting guess and apply the iteration in a loop, storing the sequence of iterates in a vector.

```{code-cell}
x = [1; zeros(4)]
for k = 1:4
    x[k+1] = x[k] - f(x[k]) / df_dx(x[k])
end
x
```

Here is the sequence of errors.

```{code-cell}
ϵ = @. x - r
```

Because the error reaches machine epsilon so rapidly, we're going to use extended precision to allow us to take a few more iterations. We'll take the last iteration as the most accurate root estimate.
```{tip}
A `BigFloat` uses 256 bits of precision, rather than 53 in `Float64`. But arithmetic is done by software emulation and is much slower.
```

```{code-cell}
x = [BigFloat(1); zeros(7)]
for k = 1:7
    x[k+1] = x[k] - f(x[k]) / df_dx(x[k])
end
r = x[end]
```

```{code-cell}
ϵ = @. Float64(x[1:end-1] - r)
```

The exponents in the scientific notation definitely suggest a squaring sequence. We can check the evolution of the ratio in {eq}`quadratictest`.

```{code-cell}
logerr = @. log10(abs(ϵ))
ratios = [NaN; [logerr[i+1] / logerr[i] for i in 1:length(logerr)-1]]
@pt :header=["iteration", "error", "log error", "ratio"] [1:7 ϵ logerr ratios]
```

The clear convergence to 2 above constitutes good evidence of quadratic convergence.
