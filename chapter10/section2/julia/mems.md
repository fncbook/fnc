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
[**Demo %s**](#demo-shooting-mems)

We revisit {numref}`Demo {number} <demo-shooting-naive>` but let {numref}`Function {number} <function-shoot>` do the heavy lifting.

```{code-cell}
λ = 0.6
ϕ = (r, w, dwdr) -> λ / w^2 - dwdr / r;
a, b = eps(), 1.0;
```

We specify the given and unknown endpoint values.

```{code-cell}
g₁(w, dw) = dw       # w' = 0 at left
g₂(w, dw) = w - 1    # w = 1 at right
r, w, dw_dx = FNC.shoot(ϕ, (a, b), g₁, g₂, [0.8, 0])
plot(r, w, title = "Shooting solution", xaxis = (L"x"), yaxis = (L"w(x)"))
```

The value of $w$ at $r=1$, meant to be exactly one, was computed to be

```{code-cell}
@show w[end];
```

The accuracy is consistent with the error tolerance used for the IVP solution. The initial value $w(0)$ that gave this solution is

```{code-cell}
@show w[1];
```
