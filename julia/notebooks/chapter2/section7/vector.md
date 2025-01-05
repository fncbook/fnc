---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---

```{index} ! Julia; norm
```

In Julia the `LinearAlgebra` package has a `norm` function for vector norms.

```{code-cell}
x = [2, -3, 1, -1]
twonorm = norm(x)         # or norm(x,2)
```

```{code-cell}
infnorm = norm(x, Inf)
```

```{code-cell}
onenorm = norm(x, 1)
```

```{index} ! Julia; normalize
```

There is also a `normalize` function that divides a vector by its norm, making it a unit vector.

```{code-cell}
normalize(x, Inf)
```
