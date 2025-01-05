---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{index} ! Julia; sprandsym
```

Here is the adjacency matrix of a graph representing a small-world network, featuring connections to neighbors and a small number of distant contacts.

```{code-cell}
using GraphRecipes
@load "smallworld.jld2" A
graphplot(A, linealpha=0.5)
```

Because each node connects to relatively few others, the adjacency matrix is quite sparse.

```{code-cell}
spy(A, title="Nonzero locations", m=2, color=:blues)
```

By {numref}`Theorem {number} <theorem-insight-adjmat>`, the entries of $\mathbf{A}^k$ give the number of walks of length $k$ between pairs of nodes, as with "*k* degrees of separation" within a social network. As $k$ grows, the density of $\mathbf{A}^k$ also grows.

```{code-cell}
plt = plot(layout=(1, 3), legend=:none, size=(600, 240))
for k in 2:4
    spy!(A^k;
        subplot=k - 1, color=:blues,
        title=latexstring("\\mathbf{A}^$k"))
end
plt
```
