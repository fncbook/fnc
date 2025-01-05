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
[**Demo %s**](#demo-qr-qrfact)


Julia provides access to both the thin and full forms of the QR factorization.

```{code-cell}
A = rand(1.:9., 6, 4)
@show m,n = size(A);
```

Here is a standard call:

```{code-cell}
Q,R = qr(A);
Q
```

```{code-cell}
R
```

If you look carefully, you see that we seemingly got a full $\mathbf{Q}$ but a thin $\mathbf{R}$. However, the $\mathbf{Q}$ above is not a standard matrix type. If you convert it to a true matrix, then it reverts to the thin form.
```{tip}
:class: dropdown
To enter the accented character `Q̂`, type `Q\hat` followed by <kbd>Tab</kbd>.
```

```{code-cell}
Q̂ = Matrix(Q)
```

We can test that $\mathbf{Q}$ is an orthogonal matrix:

```{code-cell}
opnorm(Q' * Q - I)
```

The thin $\hat{\mathbf{Q}}$ cannot be an orthogonal matrix, because it is not square, but it is still ONC:

```{code-cell}
Q̂' * Q̂ - I
```
