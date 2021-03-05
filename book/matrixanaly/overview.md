# Matrix analysis

:::{epigraph}
Judge me by my size, do you?

--Yoda, *The Empire Strikes Back* 
:::

In previous chapters we have seen how matrices that represent square or overdetermined linear systems of equations can be manipulated into LU and QR factorizations. But matrices have other factorizations that are more intrinsic to their nature as mathematical linear transformations. The most fundamental of these are the eigenvalue and singular value decompositions.

These decompositions can be used to solve linear and least squares systems, but they have greater value in how they represent the matrix itself. They lead to critical and quantitative insights about the structure of the underlying transformation and suggest ways to approximate it efficiently. In this chapter we will look at both of these fundamental decompositions and hint at just a few of their computational applications.

**Important terms**

```{glossary}
adjacency matrix
  Matrix whose nonzero entries show the links between nodes in a graph.

diagonalizable
  Matrix that admits an eigenvalue decomposition.

eigenvalue
  Scalar $\lambda$ such that $\mathbf{A}\mathbf{x} = \lambda \mathbf{x}$ for a square matrix $\mathbf{A}$ and nonzero vector $\mathbf{x}$.

eigenvalue decomposition
  Expression of a square matrix as the product of eigenvector and diagonal eigenvalue matrices.

eigenvector
  Nonzero vector $\mathbf{x}$ such that $\mathbf{A}\mathbf{x} = \lambda \mathbf{x}$ for a square matrix $\mathbf{A}$ and scalar $\lambda$.

graph
  Representation of a network as a set of nodes and edges.

hermitian
  Combination of transpose and elementwise complex conjugation. Also describes a matrix that equals its own hermitian.

hermitian positive definite
  Variant of symmetric positive definite for complex-valued matrices.

normal
  Matrix that has a unitary eigenvalue decomposition.

Rayleigh quotient
  Scalar-valued function that maps an approximate eigenvector to a result even closer to its associated eigenvalue. 

singular value decomposition
  Expression of a matrix as a product of two orthogonal/unitary matrices and a nonnegative diagonal matrix.

thin SVD
  Variant of the singular value decomposition that discards information not needed to fully represent the original matrix.

unitary
  Square matrix with complex-valued entries whose columns are orthogonal.

```

**Important Julia commands and keywords**

`eigen`
  Compute the eigenvalue decomposition of a matrix (i.e., all eigenvalues and eigenvectors).

`eigvals`
  Compute the eigenvalues of a matrix.

`graphplot`
  Plot a representation of a graph from its adjacency matrix.

`Gray`
  Convert a color value to grayscale intensity in the range $[0,1]$.

`svd`
  Compute the singular value decomposition of a matrix.

`svdvals`
  Compute the singular values of a matrix.

```{glossary}
```
