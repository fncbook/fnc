---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---

We begin with a symmetric $\mathbf{A}$.

```{code-cell}
A₁ = [  2     4     4     2
        4     5     8    -5
        4     8     6     2
        2    -5     2   -26 ];
```

We won't use pivoting, so the pivot element is at position (1,1). This will become the first element on the diagonal of $\mathbf{D}$. Then we divide by that pivot to get the first column of $\mathbf{L}$.

```{code-cell}
L = diagm(ones(4))
d = zeros(4)
d[1] = A₁[1, 1]
L[:, 1] = A₁[:, 1] / d[1]
A₂ = A₁ - d[1] * L[:, 1] * L[:, 1]'
```

We are now set up the same way for the submatrix in rows and columns 2–4.

```{code-cell}
d[2] = A₂[2, 2]
L[:, 2] = A₂[:, 2] / d[2]
A₃ = A₂ - d[2] * L[:, 2] * L[:, 2]'
```

We continue working our way down the diagonal.

```{code-cell}
d[3] = A₃[3, 3]
L[:, 3] = A₃[:, 3] / d[3]
A₄ = A₃ - d[3] * L[:, 3] * L[:, 3]'
d[4] = A₄[4, 4]
@show d;
L
```

We have arrived at the desired factorization, which we can validate:

```{code-cell}
opnorm(A₁ - (L * diagm(d) * L'))
```
