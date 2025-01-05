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
[**Demo %s**](#demo-structure-sparse)

Here we load the adjacency matrix of a graph with 2790 nodes. Each node is a web page referring to Roswell, NM, and the edges represent links between web pages. (Credit goes to Panayiotis Tsaparas and the University of Toronto for making this data public.)

```{code-cell}
using SparseArrays, JLD2
@load "roswell.jld2" A;      # file is on the book's website
```

```{index} ! Julia; nnz
```
We may define the density of $\mathbf{A}$ as the number of nonzeros divided by the total number of entries.
```{tip}
:class: dropdown
Use `nnz` to count the number of nonzeros in a sparse matrix.
```

```{code-cell}
m, n = size(A)
@show density = nnz(A) / (m * n);
```

```{index} ! Julia; summarysize
```

The computer memory consumed by any variable can be discovered using `summarysize`. We can use it to compare the space needed for the sparse representation to its dense counterpart, that is, the space needed to store all the elements, whether zero or not.

```{code-cell}
F = Matrix(A)
Base.summarysize(F) / Base.summarysize(A)
```

As you can see, the storage savings are dramatic. Matrix-vector products are also much faster using the sparse form because operations with structural zeros are skipped.

```{code-cell}
x = randn(n)
A * x;   # make sure * is loaded and compiled
@elapsed for i in 1:300
    A * x
end
```

```{code-cell}
F * x;
@elapsed for i in 1:300
    F * x
end
```
