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
[**Demo %s**](#demo-structure-sparse)

Here we load the adjacency matrix of a graph with 2790 nodes. Each node is a web page referring to Roswell, NM, and the edges represent links between web pages. (Credit goes to Panayiotis Tsaparas and the University of Toronto for making this data public.)

```{code-cell}
load roswelladj
a = whos('A')
```

```{index} ! MATLAB; nnz
```
We may define the density of $\mathbf{A}$ as the number of nonzeros divided by the total number of entries.
```{tip}
:class: dropdown
Use `nnz` to count the number of nonzeros in a sparse matrix.
```

```{code-cell}
sz = size(A);  n = sz(1);
density = nnz(A) / prod(sz)
```

The computer memory consumed by any variable can be discovered using `whos`. We can use it to compare the space needed for the sparse representation to its dense counterpart, that is, the space needed to store all the elements, whether zero or not.

```{code-cell}
F = full(A);
f = whos('F');
storage_ratio = f.bytes / a.bytes
```

Matrix-vector products are also much faster using the sparse form because operations with structural zeros are skipped.

```{code-cell}
x = randn(n,1);
tic, for i = 1:200, A*x; end
sparse_time = toc
```

```{code-cell}
tic, for i = 1:200, F*x; end
dense_time = toc
```

However, the sparse storage format in MATLAB is column-oriented.  Operations on rows may take a lot longer than similar ones on columns.

```{code-cell}
v = A(:, 1000);
tic, for i = 1:n, A(:, i) = v; end
column_time = toc
r = v';
tic, for i = 1:n, A(i, :) = r; end
row_time = toc
