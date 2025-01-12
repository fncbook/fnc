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
[**Demo %s**](#demo-structure-fill)


Here is the adjacency matrix of a graph representing a small-world network, featuring connections to neighbors and a small number of distant contacts.

```{code-cell}
load smallworld.mat
G = graph(A);
plot(G, nodecolor='r')
```

Because each node connects to relatively few others, the adjacency matrix is quite sparse.

```{code-cell}
spy(A)
```

By {numref}`Theorem {number} <theorem-insight-adjmat>`, the entries of $\mathbf{A}^k$ give the number of walks of length $k$ between pairs of nodes, as with "*k* degrees of separation" within a social network. As $k$ grows, the density of $\mathbf{A}^k$ also grows.

```{code-cell}
clf
tiledlayout(2, 2)
for k = [2, 3, 4, 6]
    nexttile
    spy(A^k)
    title(sprintf("A^{%d}", k))
end
```
