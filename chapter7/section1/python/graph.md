---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
exec(open("../../../python/FNC_init.py").read())
```
[**Demo %s**](#demo-insight-graph)

Here we create an adjacency matrix for a graph on four nodes.

```{code-cell}
A = array([[0, 1, 0, 0], [1, 0, 0, 0], [1, 1, 0, 1], [0, 1, 1, 0]])
print(A)
```

The `networkx` package has many functions for working with graphs. Here, we instruct it to create a directed graph from the adjacency matrix, then make a drawing of it.

```{index} ! Python; networkx
```

```{code-cell}
import networkx as nx
G = nx.from_numpy_array(A, create_using=nx.DiGraph)
nx.draw(G, with_labels=True, node_color="yellow")
```

Here are the counts of all walks of length 3 in the graph:

```{code-cell}
print(A**3)
```

If the adjacency matrix is symmetric, the result is an undirected graph: all edges connect in both directions.

```{code-cell}
A = array([[0, 1, 1, 0], [1, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
G = nx.from_numpy_array(A, create_using=nx.Graph)
nx.draw(G, with_labels=True, node_color="yellow")
```
