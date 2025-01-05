---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
Let's start with a known set of eigenvalues and an orthogonal eigenvector basis.

```{code-cell}
D = diagm([-6, -1, 2, 4, 5])
V, R = qr(randn(5, 5))    # V is unitary
A = V * D * V'
```

```{code-cell}
eigvals(A)
```

Now we will take the QR factorization and just reverse the factors.

```{code-cell}
Q, R = qr(A)
A = R * Q;
```

It turns out that this is a similarity transformation, so the eigenvalues are unchanged.

```{code-cell}
eigvals(A)
```

What's remarkable, and not elementary, is that if we repeat this transformation many times, the resulting matrix converges to $\mathbf{D}$.

```{code-cell}
for k in 1:40
    Q, R = qr(A)
    A = R * Q
end
A
```
