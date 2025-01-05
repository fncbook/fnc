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
[**Demo %s**](#demo-nonlinear-allencahn)


```{code-cell}
ϕ = (x, u, dudx) -> (u^3 - u) / ϵ
g₁(u, du) = du
g₂(u, du) = u - 1;
```

Finding a solution is easy at larger values of $\epsilon$.

```{code-cell}
ϵ = 0.05
init = collect(range(-1, 1, 141))
x, u₁ = FNC.bvp(ϕ, [0, 1], g₁, g₂, init)
plot(x, u₁;
    label=L"\epsilon = 0.05",  legend=:bottomright,
    xaxis=(L"x"),  yaxis=(L"u(x)"),
    title = "Allen–Cahn solution")
```

However, finding a good initialization is not trivial for smaller values of $\epsilon$. Note below that the iteration stops without converging to a solution.

```{code-cell}
ϵ = 0.002;
x, z = FNC.bvp(ϕ, [0, 1], g₁, g₂, init);
```

The iteration succeeds if we use the first solution instead as the initialization here.

```{code-cell}
x, u₂ = FNC.bvp(ϕ, [0, 1], g₁, g₂, u₁)
plot!(x, u₂; label = L"\epsilon = 0.002")
```

In this case we can continue further.

```{code-cell}
ϵ = 0.0005
x, u₃ = FNC.bvp(ϕ, [0, 1], g₁, g₂, u₂)
plot!(x, u₃, label = L"\epsilon = 0.0005")
```
