---
jupytext:
  cell_metadata_filter: -all
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: '0.8'
    jupytext_version: 1.4.2
kernelspec:
  display_name: Julia (fast start)
  language: julia
  name: julia-fast
---

# Surprising arithmetic

There is no double precision number between $1$ and $1+\epsilon_\text{mach}$. Thus the following difference is zero despite its appearance.

```{code-cell}
e = eps()/2
(1.0 + e) - 1.0
```

However, the spacing between floats in $[1.2,1)$ is $\text{mach}/2$, so both $1-\epsilon_\text{mach}/2$ and its negative are represented exactly:

```{code-cell}
1.0 + (e - 1.0)
```

This is now the "correct" result. But we have found a rather shocking breakdown of the associative law of addition!
