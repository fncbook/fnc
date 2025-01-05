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
[**Demo %s**](#demo-tensorprod-disksphere)

For a function given in polar form, such as $f(r,\theta)=1-r^4$, construction of a function over the unit disk is straightforward using a grid in $(r,\theta)$ space.

```{code-cell}
r = range(0, 1, 41)
θ = range(0, 2π, 81)
F = [1 - r^4 for r in r, θ in θ]
plot(r, θ, F';
    legend=:none, 
    color=:viridis,  fill=true,
    xlabel="r",  ylabel="θ", 
    title="A polar function")
```

Of course, we are used to seeing such plots over the $(x,y)$ plane, not the $(r,\theta)$ plane. 

In such functions the values along the line $r=0$ must be identical, and the values on the line $\theta=0$ should be identical to those on $\theta=2\pi$. Otherwise the interpretation of the domain as the unit disk is nonsensical. If the function is defined in terms of $x$ and $y$, then those can be defined in terms of $r$ and $\theta$ using {eq}`unitdiskparam`.

