---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---
```{code-cell}
:tags: [remove-cell]
include("../../../julia/FNC_init.jl")
```
[**Demo %s**](#demo-methodlines-heatBE)

Now we apply backward Euler to the heat equation. We will reuse the setup from {numref}`Demo {number} <demo-methodlines-heatFE>`. Since the matrix in {eq}`BExx` never changes during the time stepping, we do the necessary LU factorization only once.

```{code-cell}
using SparseArrays
B = sparse(I - τ * Dxx)
factor = lu(B)
for j in 1:n
    U[:, j+1] = factor \ U[:, j]
end
```

```{code-cell}
:tags: [hide-input]
idx = 1:600:n+1
times = round.(t[idx], digits=4)
label = reshape(["t = $t" for t in times], 1, length(idx))
plot(x,U[:, idx];
    label, legend=:topleft,
    title="Heat equation by backward Euler",
    xaxis=(L"x"),  yaxis=(L"u(x,0)", [0, 1]))
```

```{code-cell}
:tags: [hide-input]
anim = @animate for j in 1:20:n+1
    plot(x, U[:, j];
    label=@sprintf("t=%.5f", t[j]),
    xaxis=(L"x"),  yaxis=(L"u(x,t)", [0, 1]),
    dpi=150,  title="Heat equation by backward Euler")
end
mp4(anim, "diffusionBE.mp4")
```

This solution looks physically plausible, as the large concentration in the center diffuses outward until the solution is essentially constant. Observe that the solution remains periodic in space for all time.
