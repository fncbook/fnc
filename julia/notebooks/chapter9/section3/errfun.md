---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
We plot $|\Phi(x)|$ over the interval $[-1,1]$ with equispaced nodes for different values of $n$. 

```{code-cell} 
:tags: [hide-input]
plot(xaxis=(L"x"), yaxis=(:log10, L"|\Phi(x)|", [1e-25, 1]), legend=:bottomleft)
x = range(-1, 1, 2001)
for n in 10:10:50
    t = range(-1, 1, n+1)
    Φ(x) = prod(x - t for t in t)
    scatter!(x, abs.(Φ.(x)), m=(1, stroke(0)), label="n=$n")
end
title!("Error indicator for equispaced nodes")
```

Each time $\Phi$ passes through zero at an interpolation node, the value on the log scale should go to $-\infty$, which explains the numerous cusps on the curves.
