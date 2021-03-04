# Krylov methods in linear algebra

:::{epigraph}
I warn you not to underestimate my powers.

-— Luke Skywalker, *Return of the Jedi* 
:::

What are the implications of the $O(n^3)$ work requirements for solving linear systems? Suppose tomorrow your computer became a thousand times faster. (Historically this has taken about 15 years in the real world.) Assuming you are willing to wait just as long today as you were yesterday, the size of the linear system you can solve has gone up only by a factor of 10. Nice, but not nearly the jump that you got in hardware power. In fact, there is an odd paradox: faster computers make faster algorithms *more* important, not less, because they demand that you work at larger values of $n$, where asymptotic differences are large.

In practice the only reasonable way to deal with large matrices (at this writing, $n>10^4$ or so) is if they are sparse, or can be approximated sparsely. But LU factorization of a sparse matrix does not necessarily lead to sparse factors, particularly when row pivoting is required. The algorithm can be improved to be more sparse-aware, but we will not go into the details.

Instead, we will replace LU factorization with an iterative algorithm. Unlike the LU factorization, iteration gives useful intermediate and continually improving results before the exact solution is found, allowing us to stop well before the nominal exact termination. More importantly, though, these iterations, based on an idea called *Krylov subspaces*, allow us to fully exploit sparsity.

Krylov subspace methods have two other advantages that are subtle but critically relevant to applications. One is that they allow us to do linear algebra *even without having the relevant matrix*. This may sound undesirable or even impossible, but it exploits the [connection](matrixfree.md) between matrix-vector multiplication and a linear transformation. The other major unique advantage of Krylov subspace iterations is that they can exploit "approximate inverses" when they are available. These two features are among the most powerful ideas behind scientific computation today.

## Important terms

```{glossary}
Arnoldi iteration
  Stable algorithm for finding orthonormal bases of nested Krylov subspaces.

dominant eigenvalue
  Eigenvalue with the largest modulus (absolute value, in the real case).

fill-in
  Tendency for a sparse matrix to lose sparsity when algebraic operations are performed on it.

GMRES
  Iterative solution of a linear system through stable least-squares solutions on nested Krylov subspaces.

inverse iteration
  Use of a shift, then inverse in power iteration to transform the eigenvalue closest to a target value into a dominant one. 

Krylov matrix
  Concatenation of a vector $\mathbf{u}$ with increasing powers of a matrix times $\mathbf{u}$.

Lanczos iteration
  Specialization of the Arnoldi iteration to the case of a hermitian (or real symmetric) matrix.

power iteration
  Repeated application of a matrix to a vector, followed by normalization, resulting in convergence to an eigenvector for the dominant eigenvalue.

restarting
  Technique used in GMRES to prevent the work per iteration and overall storage from growing unboundedly.

sparse
  Term for a matrix with elements that are mostly zero for structural reasons.

upper Hessenberg
  Matrix that has nonzeros only in the upper triangle and first subdiagonal.

```

## Important Julia terms

```{glossary}
backslash
  The `\` character, which is used to solve linear systems efficiently.

```
