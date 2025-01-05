---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{index} ! Julia; enumerate
```

Suppose we want to evaluate the inverse of the function $h(x)=e^x-x$. This means solving $y=e^x-x$ for $x$ when $y$ is given, which has no elementary form. If a value of $y$ is given numerically, though, we simply have a rootfinding problem for $f(x)=e^x-x-y$.
```{tip}
The `enumerate` function produces a pair of values for each iteration: a positional index and the corresponding contents.
```

```{code-cell}
g(x) = exp(x) - x
dg_dx(x) = exp(x) - 1
y = range(g(0), g(2), 200)
x = zeros(length(y))
for (i, y) in enumerate(y)
    f(x) = g(x) - y
    df_dx(x) = dg_dx(x)
    r = FNC.newton(f, df_dx, y)
    x[i] = r[end]
end

plot(g, 0, 2, aspect_ratio=1, label=L"g(x)")
plot!(y, x, label=L"g^{-1}(y)", title="Function and its inverse")
plot!(x -> x, 0, maximum(y), label="", l=(:dash, 1), color=:black)
```
