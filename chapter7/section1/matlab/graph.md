---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
cd  /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init
```
[**Demo %s**](#demo-insight-graph)

Here we create an adjacency matrix for a graph on four nodes.

```{code-cell}
A = [0 1 0 0; 1 0 0 0; 1 1 0 1; 0 1 1 0]
```

```{index} ! MATLAB; graph (network)
```

Since this adjacency matrix is not symmetric, the edges are all directed. We use `digraph` to create a directed graph.

```{code-cell}
G = digraph(A);
plot(G)
```
Here are the counts of all walks of length 3 in the graph:

```{code-cell}
A^3
```

If the adjacency matrix is symmetric, the result is an undirected graph: all edges connect in both directions.

```{code-cell}
A = [0 1 1 0; 1 0 0 1; 1 0 0 0; 0 1 0 0];
plot(graph(A))
```

A "buckyball" is an allotrope of carbon atoms with the same connection structure as a soccer ball.

```{code-cell}
plot(graph(bucky))
```
