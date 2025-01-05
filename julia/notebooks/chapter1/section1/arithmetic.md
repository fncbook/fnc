---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---

There is no double-precision number between $1$ and $1+\epsilon_\text{mach}$. Thus the following difference is zero despite its appearance.

```{code-cell}
e = eps()/2
(1.0 + e) - 1.0
```

However, the spacing between floats in $[1/2,1)$ is $\macheps/2$, so both $1-\macheps/2$ and its negative are represented exactly:

```{code-cell}
1.0 + (e - 1.0)
```

This is now the expected result. But we have found a rather shocking breakdown of the associative law of addition!
