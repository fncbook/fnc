---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
Deleting the last row and column places all the eigenvalues of the discretization into the left half of the complex plane. 

```{code-cell}
x, Dₓ, _ = FNC.diffcheb(40, [0, 1])
A = Dₓ[1:end-1, 1:end-1];     # delete last row and column
λ = eigvals(A);
```

```{code-cell}
:tags: [hide-input]

scatter(real(λ), imag(λ);
    m=3,  aspect_ratio=1,
    legend=:none,  frame=:zerolines,
    xaxis=([-300, 100], "Re λ"),  yaxis=("Im λ"),
    title="Eigenvalues of advection with zero inflow") 
```

Note that the rightmost eigenvalues have real part at most

```{code-cell}
maximum(real(λ))
```

Consequently all solutions decay exponentially to zero as $t\to\infty$. This matches our observation of the solution: eventually, everything flows out of the domain.

