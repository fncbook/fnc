---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
We revisit {numref}`Demo %s <demo-fp-spiral>` and investigate the observed convergence more closely. Recall that above we calculated $g'(p)\approx-0.42$ at the convergent fixed point.

```{code-cell}
p = Polynomial([3.5, -4, 1])
r = roots(p)
rmin, rmax = extrema(r)
@show rmin, rmax;
```

Here is the fixed-point iteration. This time we keep track of the whole sequence of approximations.

:::{index} Julia; push!
:::

```{code-cell}
g(x) = x - p(x)
x = [2.1]
for k = 1:12
    push!(x, g(x[k]))
end
x
```

It's illuminating to construct and plot the sequence of errors.

```{code-cell}
err = @. abs(x - rmax)
plot(0:12, err;
    m=:o,
    xaxis=("iteration number"),  yaxis=("error", :log10),
    title="Convergence of fixed-point iteration")
```

It's quite clear that the convergence quickly settles into a linear rate. We could estimate this rate by doing a least-squares fit to a straight line. Keep in mind that the values for small $k$ should be left out of the computation, as they don't represent the linear trend.

```{code-cell}
y = log.(err[5:12])
p = Polynomials.fit(5:12, y, 1)
```

We can exponentiate the slope to get the convergence constant $\sigma$.

```{code-cell}
σ = exp(p.coeffs[2])
```

The error should therefore decrease by a factor of $\sigma$ at each iteration. We can check this easily from the observed data.

```{code-cell}
[err[i+1] / err[i] for i in 8:11]
```

The methods for finding $\sigma$ agree well.
