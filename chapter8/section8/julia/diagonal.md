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
[**Demo %s**](#demo-precond-diagonal)

Here is an SPD matrix that arises from solving partial differential equations.

```{code-cell}
using MatrixDepot, SparseArrays
A = matrixdepot("wathen", 60)
n = size(A, 1)
@show n, nnz(A);
```

```{index} ! Julia; DiagonalPreconditioner
```

There is an easy way to use the diagonal elements of $\mathbf{A}$, or any other vector, as a diagonal preconditioner.

```{code-cell}
using Preconditioners
b = ones(n)
M = DiagonalPreconditioner(diag(A));
```

We now compare CG with and without the preconditioner.

```{code-cell}
:tags: [hide-input]
using IterativeSolvers
plain(b) = cg(A, b, maxiter=200, reltol=1e-4, log=true)
time_plain = @elapsed x, hist1 = plain(b)
prec(b) = cg(A, b, Pl=M, maxiter=200, reltol=1e-4, log=true)
time_prec = @elapsed x, hist2 = prec(b)
@show time_plain, time_prec

rr = hist1[:resnorm]
plot(0:length(rr)-1, rr / rr[1], yscale=:log10, label="plain")
rr = hist2[:resnorm]
plot!(0:length(rr)-1, rr / rr[1], yscale=:log10, label="preconditioned")
title!("Diagonal preconditioning in CG")
```

The diagonal preconditioner cut down substantially on the number of iterations. The effect on the total time is less dramatic, but this is not a large version of the problem.

