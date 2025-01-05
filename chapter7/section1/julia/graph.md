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
[**Demo %s**](#demo-insight-graph)

Here we create an adjacency matrix for a graph on four nodes.

```{code-cell}
A = [0 1 0 0; 1 0 0 0; 1 1 0 1; 0 1 1 0]
```

```{index} ! Julia; graphplot
```

The `graphplot` function makes a visual representation of this graph.

```{code-cell}
using GraphRecipes
graphplot(A, names=1:4, markersize=0.2, arrow=6)
```

Since this adjacency matrix is not symmetric, the edges are all directed, as indicated by the arrows. Here are the counts of all walks of length 3 in the graph:

```{code-cell}
A^3
```

If the adjacency matrix is symmetric, the result is an undirected graph: all edges connect in both directions.

```{code-cell}
A = [0 1 1 0; 1 0 0 1; 1 0 0 0; 0 1 0 0]
graphplot(A, names=1:4, markersize=0.2)
```
