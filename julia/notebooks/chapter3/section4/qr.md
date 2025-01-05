---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---

We will use Householder reflections to produce a QR factorization of a random matrix.
```{tip}
The `rand` function can select randomly from within the interval $[0,1]$, or from a vector or range that you specify.
```

```{code-cell}
A = rand(float(1:9), 6, 4)
m,n = size(A)
```

```{index} Julia; normalize, ! Julia; I
```

Our first step is to introduce zeros below the diagonal in column 1 by using {eq}`hhvector` and {eq}`hhreflect`. 
:::
```{tip}
`I` can stand for an identity matrix of any size, inferred from the context when needed.
```

```{code-cell}
z = A[:, 1];
v = normalize(z - norm(z) * [1; zeros(m-1)])
P₁ = I - 2v * v'   # reflector
```

We check that this reflector introduces zeros as it should:

```{code-cell}
P₁ * z
```

Now we replace $\mathbf{A}$ by $\mathbf{P}\mathbf{A}$.

```{code-cell}
A = P₁ * A
```

We are set to put zeros into column 2. We must not use row 1 in any way, lest it destroy the zeros we just introduced. So we leave it out of the next reflector.

```{code-cell}
z = A[2:m, 2]
v = normalize(z - norm(z) * [1; zeros(m-2)])
P₂ = I - 2v * v'
```

We now apply this reflector to rows 2 and below only.

```{code-cell}
A[2:m, :] = P₂ * A[2:m, :]
A
```

We need to iterate the process for the last two columns.

```{code-cell}
for j in 3:n
    z = A[j:m, j]
    v = normalize(z - norm(z) * [1; zeros(m-j)])
    P = I - 2v * v'
    A[j:m, :] = P * A[j:m, :]
end
```

We have now reduced the original to an upper triangular matrix using four orthogonal Householder reflections:

```{code-cell}
R = triu(A)
```
