---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering: 
  headings: false
---
Let's define a set of four nodes (i.e., $n=3$ in our formulas).

```{index} ! Julia; annotate!
```

```{code-cell}
t = [0, 0.55, 0.7, 1]
```

We plot the hat functions $H_0,\ldots,H_3$.
```{tip}
Use `annotate!` to add text to a plot.
```

```{code-cell}
using Plots
plt = plot(layout=(4, 1),  legend=:top,
    xlabel=L"x",  ylims=[-0.1, 1.1],  ytick=[])
for k in 0:3
    Hₖ = FNC.hatfun(t, k)
    plot!(Hₖ, 0, 1, subplot=k + 1)
    scatter!(t, Hₖ.(t), m=3, subplot=k + 1)
    annotate!(t[k+1], 0.25, text(latexstring("H_$k"), 10), subplot=k+1)
end
plt
```
