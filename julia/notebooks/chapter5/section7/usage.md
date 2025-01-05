---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering: 
  headings: false
---
We'll integrate the function from {numref}`Demo %s <demo-adapt-motive>`.

```{code-cell}
f = x -> (x + 1)^2 * cos((2x + 1) / (x - 4.3));
```

We perform the integration and show the nodes selected underneath the curve.

```{code-cell}
A, t = FNC.intadapt(f, 0, 4, 0.001)
@show num_nodes = length(t);

plot(f, 0, 4;
    color=:black, legend=:none,
    xlabel=L"x",  ylabel=L"f(x)", 
    title="Adaptive node selection")
plot!(t, f.(t), seriestype=:sticks, m=(:o, 2))
```

The error turns out to be a bit more than we requested. It's only an estimate, not a guarantee.

```{code-cell}
Q, _ = quadgk(f, 0, 4, atol=1e-14, rtol=1e-14);    # 'exact' value
println("error: $(Q-A)");
```

Let's see how the number of integrand evaluations and the error vary with the requested tolerance.

```{code-cell}
tol = [1 / 10^k for k in 4:14]
err, n = [], []
for tol in 10.0 .^ (-4:-1:-14)
    A, t = FNC.intadapt(f, 0, 4, tol)
    push!(err, Q - A)
    push!(n, length(t))
end
@pt :header=["tolerance", "error", "number of nodes"] [tol err n][1:2:end, :]
```

As you can see, even though the errors are not smaller than the tolerances, the two columns decrease in tandem. If we consider now the convergence not in $h$, which is poorly defined now, but in the number of nodes actually chosen, we come close to the fourth-order accuracy of the underlying Simpson scheme.

```{code-cell}
plot(n, abs.(err);
    m=:o, label="results",
    xaxis=(:log10, "number of nodes"),  yaxis=(:log10, "error"),
    title="Convergence of adaptive integration")

order4 = @. 0.01 * (n / n[1])^(-4)
plot!(n, order4, l=:dash, label=L"O(n^{-4})")
```
