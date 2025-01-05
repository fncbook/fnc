---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{code-cell}
include("../../../julia/FNC_init.jl")
```
[**Demo %s**](#demo-precond-gmres)

Here is a nonsymmetric matrix arising from a probabilistic model in computational chemistry.

```{code-cell}
using SparseArrays, MatrixDepot
A = sparse(matrixdepot("Watson/chem_master1"))
n = size(A, 1)
@show n, nnz(A), issymmetric(A)
```

Without a preconditioner, GMRES makes essentially no progress after 100 iterations.

```{code-cell}
b = rand(40000)
using IterativeSolvers
const GMRES = IterativeSolvers.gmres
x, history = GMRES(A, b, maxiter=100, reltol=1e-5, log=true)
resnorm = history[:resnorm]
@show resnorm[end] / resnorm[1];
```

```{index} ! Julia; ilu
```

The following version of incomplete LU factorization drops all sufficiently small values (i.e., replaces them with zeros). This keeps the number of nonzeros in the factors under control.

```{code-cell}
using IncompleteLU
iLU = ilu(A, τ=0.25)
@show nnz(iLU) / nnz(A);
```

The result is almost 10 times as dense as $\mathbf{A}$ and yet still not a true factorization of it. However, it's close enough for an approximate inverse in a preconditioner. The actual preconditioning matrix is $\mathbf{M}=\mathbf{L}\mathbf{U}$, but we just supply the factorization to `gmres`.

```{code-cell}
_, history = GMRES(A, b, Pl=iLU, maxiter=100, reltol=1e-5, log=true)
history
```

The $\tau$ parameter in `ilu` balances the accuracy of the iLU factorization with the time needed to compute it and invert it. As $\tau\to 0$, more of the elements are kept, making the preconditioner more effective but slower per iteration.

```{code-cell}
:tags: [hide-input]

plt = plot(0:40, resnorm[1:41] / resnorm[1];
    label="no preconditioning",  legend=:bottomright,
    xaxis=("iteration number"),
    yaxis=(:log10, "residual norm"),
    title="Incomplete LU preconditioning")
for τ in [2, 1, 0.25, 0.1]
    t = @elapsed iLU = ilu(A; τ)
    t += @elapsed _, history = GMRES(A, b, Pl=iLU, maxiter=100,
        reltol=1e-5, log=true)
    resnorm = history[:resnorm]
    label = "τ = $τ, time = $(round(t,digits=3))"
    plot!(0:length(resnorm)-1, resnorm / resnorm[1]; label)
end
plt
```

In any given problem, it's impossible to know in advance where the right balance lies between fidelity and speed for the preconditioner.
