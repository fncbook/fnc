---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
include("../../../julia/FNC_init.jl")
```
[**Demo %s**](#demo-linear-solve)

```{code-cell}
exact = x -> exp(sin(x));
```

The problem is presented above in our standard form, so we can identify the coefficient functions in the ODE. Each should be coded as a function.

```{code-cell}
p = x -> -cos(x);
q = sin;
r = x -> 0;      # function, not value 
```

We solve the BVP and compare the result to the exact solution.

```{code-cell}
x, u = FNC.bvplin(p, q, r, [0, 3π / 2], 1, exp(-1), 30);
```

```{code-cell}
plot(exact, 0, 3π / 2, layout = (2, 1), label = "exact")
scatter!(x, u, m = :o,
    subplot=1,  label="numerical",
    yaxis=("solution"),
    title="Solution of a linear BVP")
plot!(x, exact.(x) - u, subplot = 2, xaxis = L"x", yaxis = ("error"))
```
