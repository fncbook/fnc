---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
Here are the data for a linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$. 

```{code-cell}
A = [2 0 4 3; -4 5 -7 -10; 1 15 2 -4.5; -2 0 2 -13];
b = [4,9,9,4];
```

We apply {numref}`Function {number} <function-lufact>` and then do two triangular solves.

```{code-cell}
L, U = FNC.lufact(A)
z = FNC.forwardsub(L, b)
x = FNC.backsub(U, z)
```

A check on the residual assures us that we found the solution.

```{code-cell}
b - A*x
```
