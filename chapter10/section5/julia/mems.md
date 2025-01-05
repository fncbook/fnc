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
[**Demo %s**](#demo-nonlinear-mems)

Here is the problem definition. We use a truncated domain to avoid division by zero at $r=0$.

```{code-cell}
domain = [eps(), 1]
λ = 0.5
ϕ = (r, w, dwdr) -> λ / w^2 - dwdr / r
g₁(w, dw) = dw
g₂(w, dw) = w - 1;
```

First we try a constant function as the initialization.

```{code-cell}
init = ones(301)
r, w₁ = FNC.bvp(ϕ, domain, g₁, g₂, init)

plot(r, w₁;
    xaxis = (L"r"),  yaxis = (L"w(r)"), 
    title = "Solution of the MEMS problem")
```

It's not necessary that the initialization satisfy the boundary conditions. In fact, by choosing a different constant function as the initial guess, we arrive at another valid solution.

```{code-cell}
init = 0.5 * ones(301)
r, w₂ = FNC.bvp(ϕ, domain, g₁, g₂, init)
plot!(r, w₂, title = "Two solutions of the MEMS problem")
```
