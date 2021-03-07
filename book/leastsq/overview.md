# Overdetermined linear systems

```{epigraph}
I must have hit her pretty close to the mark to get her all riled up like that, huh, kid? 

-- Han Solo, *The Empire Strikes Back*
```

So far we have considered $\mathbf{A}\mathbf{x}=\mathbf{b}$ only when $\mathbf{A}$ is a square matrix. In this chapter we consider how to interpret and solve the problem for an $m\times n$ matrix where $m>n$—and in practice, $m$ is often *much* larger than $n$. This is called an {term}`overdetermined` linear system because, in general, the system has more equations to satisfy than the variables allow. The complementary *underdetermined* case $m<n$ turns up less frequently and will not be considered in this book.

Since we cannot solve all of the system's equations, we need to define what the "best possible" answer is. There are multiple useful options, but the most important version of the overdetermined problem occurs using the **least squares**—the sum of the squares of the equation residuals is minimized. This is far from an arbitrary choice. Mathematically, we recognize the sum-of-squares as a vector 2-norm and therefore tied to inner products; physically, the 2-norm may coincide with energy, which is often minimized by natural systems; and statistically, least squares leads to the estimates of maximum likelihood for certain models. Furthermore the solution of the least-squares problem requires only linear algebra and is about as easily to compute as in the square case.

The linear least squares problem serves as our introduction to the vast field of *optimization*. It is one of the simplest problems of this type. We will see an extension to a nonlinear version in the next chapter.

**Important terms**

```{glossary}
linear least squares problem
  Minimization of the 2-norm of the residual for an overdetermined linear system.

normal equations
  Square linear system equivalent to the linear least squares problem.

ONC matrix
  Matrix that has orthonormal columns.

orthogonal
  Vectors that have an inner product of zero.

orthogonal matrix
  Square ONC matrix, i.e., matrix whose transpose is its inverse.

orthonormal
  Vectors that are both mutually orthogonal and all of unit 2-norm.

overdetermined
  Characterized by having more constraints than available degrees of freedom.

pseudoinverse
  Rectangular matrix that maps data to solution in the linear least squares problem.

QR factorization
  Representation of a matrix as the product of an orthogonal and an upper triangular matrix.
```

**Important Julia terms**

```{glossary}
`pinv`
  Pseudoinverse of a rectangular matrix.

`qr`
  Full or thin QR factorization of a matrix.
```
