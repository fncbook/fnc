---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---

```{index} Julia; tril, Julia; triu
```
We explore the outer product formula for two random triangular matrices.

```{code-cell}
L = tril( rand(1:9, 3, 3) )
```

```{code-cell}
U = triu( rand(1:9, 3, 3) )
```

Here are the three outer products in the sum in {eq}`matrixouter`:
```{tip}
Although `U[1,:]` is a row of `U`, it is a vector, and as such it has a default column interpretation.
```

```{code-cell}
L[:, 1] * U[1, :]'
```

```{code-cell}
L[:, 2] * U[2, :]'
```

```{code-cell}
L[:, 3] * U[3, :]'
```

Simply because of the triangular zero structures, only the first outer product contributes to the first row and first column of the entire product. 
