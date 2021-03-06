# Sparsity and structure

```{index} matrix; sparse
``` 
```{index} sparse matrix
```

Very large matrices cannot be stored all within primary memory of a computer unless they are {term}`sparse`. A sparse matrix has *structural zeros*, meaning entries that are known to be exactly zero.} For instance, the adjacency matrix of a graph has zeros where there are no links in the graph. To store and operate with a sparse matrix efficiently, it is not represented as an array of all of its values. There is a variety of sparse formats available; MATLAB chooses *compressed sparse column* format. For the most part, you can imagine that the matrix is stored as triples $(i,j,A_{ij})$ for all the nonzero $(i,j)$ locations.

## Computing with sparse matrices

```{index} adjacency matrix
```

Most graphs with real applications have many fewer edges than the maximum possible $n^2$ for $n$ nodes. Accordingly, their adjacency matrices have mostly zero elements and should be represented sparsely. 

::::{prf:example} Julia demo
:class: demo
:label: demos-structure-sparse
{doc}`demos/structure-sparse`
::::

Arithmetic operations such as "+", "-", "*", and \verb?^? respect and exploit sparsity, if the matrix operands are sparse. However, matrix operations may substantially decrease the amount of sparsity, a phenomenon known as {term}`fill-in`.

In the case of an adjacency matrix $\mathbf{A}$, for example, the $(i,j)$ entry of matrix $\mathbf{A}^k$ for positive integer $k$ is the number of paths of length $k$ from node $i$ to node $j$.


::::{prf:example} Julia demo
:class: demo
:label: demos-structure-fill
{doc}`demos/structure-fill`
::::

## Banded matrices

```{index} matrix; banded
```

A particularly important type of sparse matrix is a banded matrix. [Earlier](../linsys/structure.md) we said that $\mathbf{A}$ has **upper bandwith** $p$ if $j-i > p$ implies $A_{ij}=0$, and **lower bandwidth** $q$ if $i-j > q$ implies $A_{ij}=0$. We say the total **bandwidth** is $p+q+1$. 

Banded matrices appear naturally in many applications where each node interacts directly with only its closest neighbors. Without pivoting, an LU factorization preserves bandwidth, but pivoting can change or destroy bandedness.

```{index} matrix; factorization
```

::::{prf:example} Julia demo
:class: demo
:label: demos-structure-banded
{doc}`demos/structure-banded`
::::

## Linear systems and eigenvalues

If given a sparse matrix, the backslash operator will automatically try a form of sparse-aware Cholesky or pivoted LU factorization. Depending on the sparsity pattern of the matrix, the time taken to solve the linear system may be well below the $O(n^3)$ needed in the general case.

```{index} eigenvalue decomposition
```

For very large matrices, it's unlikely that you will want to find all of its eigenvalues and eigenvectors. In an [upcoming section](subspace.md) we will describe some of the math behind an algorithm that can find a selected number of eigenvalues of largest magnitude, lying to the extreme left or right, or nearest a given complex number. 

::::{prf:example} Julia demo
:class: demo
:label: demos-structure-linalg
{doc}`demos/structure-linalg`
::::

<!-- \begin{exercises}
    \input{krylov/exercises/Sparse}
\end{exercises} -->
