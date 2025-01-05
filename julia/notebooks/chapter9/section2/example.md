---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{code-cell}
using Plots
f(x) = sin(exp(2x))
plot(f, 0, 1, label="function", legend=:bottomleft)
```

```{code-cell}
t = (0:3) / 3
y = f.(t)
scatter!(t, y, color=:black, label="nodes")
```

```{code-cell}
p = FNC.polyinterp(t, y)
plot!(p, 0, 1, label="interpolant", title="Interpolation on 4 nodes")
```

The curves must intersect at the interpolation nodes. For $n=6$ the interpolant is noticeably better.

```{code-cell}
plot(f, 0, 1, label="function", legend=:bottomleft)
t = (0:6) / 6
y = f.(t)
p = FNC.polyinterp(t, y)
scatter!(t, y, color=:black, label="nodes")
plot!(p, 0, 1, label="interpolant", title="Interpolation on 7 nodes")
```
