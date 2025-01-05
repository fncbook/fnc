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
[**Demo %s**](#demo-absstab-advdiff)

The eigenvalues of advection-diffusion are near-imaginary for $\epsilon\approx 0$ and get closer to the negative real axis as $\epsilon$ increases.

```{code-cell}
:tags: [hide-input]
plt = plot(
    legend=:topleft,
    aspect_ratio=1,
    xlabel="Re ζ",  ylabel="Im ζ",
    title="Eigenvalues for advection-diffusion")
x, Dₓ, Dₓₓ = FNC.diffper(40, [0, 1]);
for ϵ in [0.001, 0.01, 0.05]
    λ = eigvals(-Dₓ + ϵ*Dₓₓ)
    scatter!(real(λ), imag(λ), m=:o, label="\\epsilon = $ϵ")
end
plt
```
