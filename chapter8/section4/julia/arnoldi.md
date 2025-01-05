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
[**Demo %s**](#demo-subspace-arnoldi)

We illustrate a few steps of the Arnoldi iteration for a small matrix.

```{code-cell}
A = rand(1.0:9.0, 6, 6)
```

The seed vector we choose here determines the first member of the orthonormal basis.

```{code-cell}
u = randn(6)
Q = u / norm(u);
```

Multiplication by $\mathbf{A}$ gives us a new vector in $\mathcal{K}_2$.

```{code-cell}
Aq = A * Q[:, 1];
```

We subtract off its projection in the previous direction. The remainder is rescaled to give us the next orthonormal column.

```{code-cell}
v = Aq - dot(Q[:, 1], Aq) * Q[:, 1]
Q = [Q v / norm(v)];
```

On the next pass, we have to subtract off the projections in two previous directions.

```{code-cell}
Aq = A * Q[:, 2]
v = Aq - dot(Q[:, 1], Aq) * Q[:, 1] - dot(Q[:, 2], Aq) * Q[:, 2]
Q = [Q v / norm(v)];
```

At every step, $\mathbf{Q}_m$ is an ONC matrix.

```{code-cell}
@show opnorm(Q' * Q - I);
```

And $\mathbf{Q}_m$ spans the same space as the three-dimensional Krylov matrix.

```{code-cell}
K = [u A * u A * A * u];
@show rank([Q K]);
```
