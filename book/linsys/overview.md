# Linear systems of equations

```{epigraph}
It's all a lot of simple tricks and nonsense.

-- Han Solo, *Star Wars: A New Hope*
```

One of the most frequently encountered tasks in scientific computation is the solution of the linear system of equations $\mathbf{A} \mathbf{x}=\mathbf{b}$ for a given square matrix $\mathbf{A}$ and vector $\mathbf{b}$.  This problem can be solved in a finite number of steps, using an algorithm equivalent to Gaussian elimination. Describing the algorithm is mostly an exercise in organizing some linear algebra.

Analyzing the algorithm requires new tools. Because the computations will take place in floating point, we must first discuss a system for measuring the "size" of a perturbation to a vector or matrix data. Once that is understood, we find that the conditioning of the square linear system problem is quite straightforward to describe. Finally, we will see that the algorithm may change when certain things are known about the matrix $\mathbf{A}$.

**Important terms**

```{glossary}
asymptotic
  Asymptotic relationship indicating "same growth up to less significant terms" of a standard function.

backward substitution
  Systematic method for solving a linear system with an upper triangular matrix.

bandwidth
  The number of diagonals around the main diagonal that have nonzero elements.

big-O
  Asymptotic relationship indicating "growth no faster than" a standard function.

Cholesky factorization
  Symmetrized version of LU factorization for SPD matrices.

flops
  Arithmetic operations on floating-point numbers.

forward substitution
  Systematic method for solving a linear system with a lower triangular matrix.

Frobenius norm
  Matrix norm computed by applying the vector 2-norm to a linearized interpretation.

Gaussian elimination
  Use of row operations to transform a linear system to an equivalent one in triangular form.

hermitian
  Conjugate transpose of a complex matrix.

identity matrix
  Matrix with ones on the diagonal, acting as the multiplicative identity.

induced matrix norm
  Norm computed using the interpretation of a matrix as a linear operator.

interpolation
  Constructing a function that passes through a given set of data points.

LU factorization
  Factorization of a square matrix into the product of a unit lower triangular matrix and an upper triangular matrix.

matrix condition number
  Norm of the matrix times the norm of its inverse.

norm
  Means of defining the magnitude of a vector or matrix.

permutation matrix
  Identity matrix with reordered rows.

PLU factorization
  LU factorization with row pivoting.

residual
  For a linear system, the difference between $\mathbf{b}$ and $\mathbf{A}\tilde{\mathbf{x}}$ for a computed solution approximation $\tilde{\mathbf{x}}$.

row pivoting
  Swapping rows during PLU factorization to ensure that the factorization exists and can be computed stably.

sparse
  Matrix with mainly zero entries for structural reasons.

symmetric matrix
  Square matrix that is equal to its transpose.

symmetric positive definite matrix
  Matrix that is symmetric and positive definite, thereby permitting a Cholesky factorization.

triangular matrix
  Matrix that is all zero either above (for lower triangular) or below (for upper triangular) the main diagonal.

tridiagonal matrix
  Matrix with nonzeros only on the main diagonal and the adjacent two diagonals.

unit triangular matrix
  Triangular matrix that has a 1 in every position on the main diagonal.

unit vector
  A vector whose norm equals one.

Vandermonde matrix
  Type of matrix important to posing polynomial interpolation conditions.
```

**Important Julia terms**

```{glossary}
backslash
  The `\` character, which is used to solve linear systems efficiently.

comprehension
  Shorthand alternative to loops when creating some vectors and matrices.

`copy`
  Copy the contents of an array into a new array, so that changes to the new array don't affect the old one.

`end`
  Last element of an array in a particular dimension.

`float`
  Convert scalar or array to floating-point representation.

`InexactError`
  Thrown when an implied type conversion is impossible for the given value.

`inv`
  Inverse of a matrix, *not* to be used for solving linear systems.

`norm`
  Norm of a vector, or Frobenius norm of a matrix.

`ones`
  Create a vector or matrix of ones.

`opnorm`
  Matrix norm in the induced (operator) sense.

`size`
  Return the size (sometimes called shape) of a vector or matrix.

`spy`
  Visualize the nonzero structure of a sparse matrix.

`zeros`
  Create a vector or matrix of zeros.
```
