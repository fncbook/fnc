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

# Stable alternative to the quadratic formula

We repeat the rootfinding experiment of {doc}`stability-quadbad` with an alternative algorithm.

```{code-cell}
a = 1;  b = -(1e6+1e-6);  c = 1;
```

First, we find the "good" root using the quadratic formula.

```{code-cell}
@show x1 = (-b + sqrt(b^2-4*a*c)) / (2*a);
```

Then we use the alternative formula for computing the smaller root:

```{code-cell}
@show x2 = c/(a*x1);
```

As you see in this output, Julia often suppresses trailing zeros in a decimal expansion. To be sure we have an accurate result, we compute its relative error.

```{code-cell}
abs(x2-1e-6) / 1e-6
```
