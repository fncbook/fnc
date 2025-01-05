---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
We repeat the rootfinding experiment of {numref}`Demo %s <demo-stability-quadbad>` with an alternative algorithm.

```{code-cell}
a = 1;  b = -(1e6 + 1e-6);  c = 1;
```

First, we find the "good" root using the quadratic formula.

```{code-cell}
@show x₁ = (-b + sqrt(b^2 - 4a*c)) / 2a;
```

Then we use the identity $x_1x_2=\frac{c}{a}$ to compute the smaller root:

```{code-cell}
@show x₂ = c / (a * x₁);
```

As you see in this output, Julia often suppresses trailing zeros in a decimal expansion. To be sure we have an accurate result, we compute its relative error.

```{code-cell}
abs(x₂ - 1e-6) / 1e-6
```
