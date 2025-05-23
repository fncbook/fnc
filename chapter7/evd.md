---
numbering:
  enumerator: 7.2.%s
---
(section-matrixanaly-evd)=

# Eigenvalue decomposition

To this point we have dealt frequently with the solution of the linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$. Alongside this problem in its importance to linear algebra is the eigenvalue problem.

```{index} ! eigenvalue, ! eigenvector
```

::::{prf:definition} Eigenvalue and eigenvector
:label: definition-eigenvalue
Given a square matrix $\mathbf{A}$, if

```{math}
:label: eigdef
\mathbf{A}\mathbf{x} = \lambda \mathbf{x}
```

for a scalar $\lambda$ and a nonzero vector $\mathbf{x}$, then $\lambda$ is an {term}`eigenvalue` and $\mathbf{x}$ is an associated {term}`eigenvector`.
::::

## Complex matrices

A matrix with real entries can have complex eigenvalues. Therefore, we assume all matrices, vectors, and scalars may be complex in what follows. Recall that a complex number can be represented as $a+i b$ for real $a$ and $b$ and where $i^2=-1$. The **complex conjugate** of $x=a+i b$ is denoted $\bar{x}$ and is given by $\bar{x}=a-i b$. The **magnitude** or *modulus* of a complex number $z$ is

$$
|z| = \sqrt{z\cdot \bar{z}}.
$$

```{index} ! unitary matrix, orthogonal matrix
```

::::{prf:definition} Terms for complex matrices
:label: definition-complexmatrix
The {term}`adjoint` or *hermitian* of a matrix $\mathbf{A}$ is denoted $\mathbf{A}^*$ and is given by $\mathbf{A}^*=(\overline{\mathbf{A}})^T=\overline{\mathbf{A}^T}$. The matrix is **self-adjoint** or {term}`hermitian` if $\mathbf{A}^*=\mathbf{A}$.

The **2-norm** of a complex vector $\mathbf{u}$ is $\sqrt{\mathbf{u}^*\mathbf{u}}$. Other vector norms, and all matrix norms, are as defined in {numref}`section-linsys-norms`.

Complex vectors $\mathbf{u}$ and $\mathbf{v}$ of the same dimension are {term}`orthogonal vectors` if $\mathbf{u}^*\mathbf{v}=0$ and are {term}`Orthonormal vectors` if both also have unit 2-norm. A {term}`unitary matrix` is a square matrix with orthonormal columns, or, equivalently, a matrix satisfying $\mathbf{A}^* = \mathbf{A}^{-1}$.
::::

For the most part, "adjoint" replaces "transpose," "hermitian" replaces "symmetric," and "unitary matrix" replaces "orthogonal matrix" when applying our previous results to complex matrices.

## Eigenvalue decomposition

```{index} ! characteristic polynomial
```

An easy rewrite of the eigenvalue definition {eq}`eigdef` is that $(\mathbf{A} - \lambda\mathbf{I}) \mathbf{x} = \boldsymbol{0}$. Hence, $(\mathbf{A} - \lambda\mathbf{I})$ is singular, and it therefore must have a zero determinant. This is the property most often used to compute eigenvalues by hand.

::::{prf:example}
Given

```{math}
\mathbf{A} = \begin{bmatrix} 1 & 1 \\ 4 & 1 \end{bmatrix},
```

we compute

```{math}
\begin{vmatrix}
1-\lambda & 1\\ 
4 & 1-\lambda
\end{vmatrix}
= (1-\lambda)^2 - 4 = \lambda^2-2\lambda-3.
```

The eigenvalues are the roots of this quadratic, $\lambda_1=3$ and $\lambda_2=-1$.
::::

The determinant $\det(\mathbf{A} - \lambda \mathbf{I})$ is called the **characteristic polynomial**. Its roots are the eigenvalues, so we know that an $n\times n$ matrix has $n$ eigenvalues, counting algebraic multiplicity.

Suppose that $\mathbf{A}\mathbf{v}_k=\lambda_k\mathbf{v}_k$ for $k=1,\ldots,n$. We can summarize these as

```{math}
\begin{split}
   \begin{bmatrix}
    \mathbf{A}\mathbf{v}_1 & \mathbf{A}\mathbf{v}_2 & \cdots & \mathbf{A}\mathbf{v}_n
  \end{bmatrix}
  &=
    \begin{bmatrix}
      \lambda_1 \mathbf{v}_1 & \lambda_2\mathbf{v}_2 & \cdots & \lambda_n \mathbf{v}_n
    \end{bmatrix}, \\[1mm]
  \mathbf{A} \begin{bmatrix}
    \mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n
  \end{bmatrix}
  &=
\begin{bmatrix}
    \mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n
  \end{bmatrix}
  \begin{bmatrix}
    \lambda_1 & &  &  \\
    & \lambda_2 & & \\
    & & \ddots & \\
    & & & \lambda_n
  \end{bmatrix},
\end{split}
```

which we write as

```{math}
:label: ev-all
  \mathbf{A} \mathbf{V} = \mathbf{V} \mathbf{D}.
```

If we find that $\mathbf{V}$ is a nonsingular matrix, then we arrive at a key factorization.[^factdecomp]

[^factdecomp]: The terms "factorization" and "decomposition" are equivalent; they coexist mainly for historical reasons.

```{index} ! matrix factorization; EVD, diagonal matrix
```

```{index} ! eigenvalue decomposition, ! diagonalizable matrix
```

::::{prf:definition} Eigenvalue decomposition (EVD)
:label: definition-evd
An **eigenvalue decomposition** (EVD) of a square matrix $\mathbf{A}$ is

```{math}
:label: evdecomp
\mathbf{A} = \mathbf{V} \mathbf{D} \mathbf{V}^{-1}.
```

If $\mathbf{A}$ has an EVD, we say that $\mathbf{A}$ is a {term}`diagonalizable matrix`; otherwise $\mathbf{A}$ is **nondiagonalizable** (or *defective*).
::::

Observe that if $\mathbf{A}\mathbf{v} = \lambda \mathbf{v}$ for nonzero $\mathbf{v}$, then the equation remains true for any nonzero multiple of $\mathbf{v}$. Therefore:

:::{important}
A matrix that has an eigenvalue decomposition has infinitely many of them.
:::

We stress that while {eq}`ev-all` is possible for all square matrices, {eq}`evdecomp` is not.  One simple example of a nondiagonalizable matrix is

```{math}
:label: jordanblock
  \mathbf{B} = \begin{bmatrix}
    1 & 1\\0 & 1
  \end{bmatrix}.
```

There is a common circumstance in which we can guarantee an EVD exists. The proof of the following theorem can be found in many elementary texts on linear algebra.

````{prf:theorem}
If the $n\times n$ matrix $\mathbf{A}$ has $n$ distinct eigenvalues, then $\mathbf{A}$ is diagonalizable.
````

::::{prf:example}  Eigenvalues and eigenvectors
:label: demo-evd-eigen

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-evd-eigen-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-evd-eigen-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-evd-eigen-python
:::
```` 
`````

::::

## Similarity and matrix powers

```{index} ! similarity transformation
```

The particular relationship between matrices $\mathbf{A}$ and $\mathbf{D}$ in {eq}`evdecomp` is important.

:::{prf:definition} Similar matrices
:label: definition-similar
If $\mathbf{S}$ is any nonsingular matrix, we say that $\mathbf{B}=\mathbf{S}\mathbf{A}\mathbf{S}^{-1}$ is a {term}`similarity transformation` of $\mathbf{A}$, and we say that $\mathbf{A}$ and $\mathbf{B}$ are {term}`similar matrices`.
:::

Hence, an EVD transforms $\mathbf{A}$ to a similar matrix that happens to be diagonal, which is as simple as a matrix gets.

One way to interpret similarity is via change of basis (see @obs-basis):

```{math}
\mathbf{B}\mathbf{x} = \mathbf{S}\mathbf{A}\mathbf{S}^{-1} \mathbf{x} 
= \underbrace{\mathbf{S} \underbrace{ \Bigl(\mathbf{A} \underbrace{\left( \mathbf{S}^{-1} \mathbf{x}\right)}_{\text{into $S$-basis}}\Bigr)}_{\text{apply $\mathbf{A}$}}}_{\text{out of $S$-basis}} .
```

That is, $\mathbf{A}$ and $\mathbf{B}$ represent the same linear transformation in different bases.

A similarity transformation does not change eigenvalues, a fact that is typically proved in elementary linear algebra texts:

````{prf:theorem}
If $\mathbf{S}$ is a nonsingular matrix, then $\mathbf{S}\mathbf{A}\mathbf{S}^{-1}$ has the same eigenvalues as $\mathbf{A}$.
````

The EVD is especially useful for matrix powers. To begin,

```{math}
\mathbf{A}^2=(\mathbf{V}\mathbf{D}\mathbf{V}^{-1})(\mathbf{V}\mathbf{D}\mathbf{V}^{-1})=\mathbf{V}\mathbf{D}(\mathbf{V}^{-1}\mathbf{V})\mathbf{D}\mathbf{V}^{-1}=\mathbf{V}\mathbf{D}^2\mathbf{V}^{-1}.
```

Multiplying this result by $\mathbf{A}$ repeatedly, we find that

```{math}
:label: evdpower
\mathbf{A}^k = \mathbf{V}\mathbf{D}^k\mathbf{V}^{-1}.
```

Because $\mathbf{D}$ is diagonal, its power $\mathbf{D}^k$ is just the diagonal matrix of the $k$th powers of the eigenvalues.

Furthermore, given a polynomial $p(z)=c_0+c_1 z + \cdots + c_m z^m$, we can apply the polynomial to the matrix in a straightforward way,

```{math}
:label: matrixpoly
p(\mathbf{A}) = c_0\mathbf{I}  +c_1 \mathbf{A} + \cdots + c_m \mathbf{A}^m.
```

Applying {eq}`evdpower` leads to

```{math}
:label: matrixpolyevd
\begin{split}
p(\mathbf{A}) & = c_0\mathbf{V}\mathbf{V}^{-1}  +c_1 \mathbf{V}\mathbf{D}\mathbf{V}^{-1} + \cdots + c_m \mathbf{V}\mathbf{D}^m\mathbf{V}^{-1} \\ 
&= \mathbf{V} \cdot [ c_0\mathbf{I}  +c_1 \mathbf{D} + \cdots + c_m \mathbf{D}^m] \cdot \mathbf{V}^{-1} \\[1mm] 
&= \mathbf{V} \cdot \begin{bmatrix}
  p(\lambda_1) & & & \\ & p(\lambda_2) & &  \\ & & \ddots & \\ & & & p(\lambda_n)  
\end{bmatrix} \cdot \mathbf{V}^{-1}.
\end{split}
```

Finally, given the convergence of Taylor polynomials to common functions, we are able to apply a function $f$ to a square matrix by replacing $p$ with $f$ in {eq}`matrixpolyevd`.

## Conditioning of eigenvalues

```{index} ! eigenvalue; conditioning of, condition number; of eigenvalues
```

```{index} ! Bauer–Fike theorem
```

Just as linear systems have condition numbers that quantify the effect of finite precision, eigenvalue problems may be poorly conditioned too. While many possible results can be derived, we will use just one, the **Bauer–Fike theorem**.

````{prf:theorem} Bauer–Fike
:label: theorem-bauer-fike
Let $\mathbf{A}\in\mathbb{C}^{n\times n}$ be diagonalizable, $\mathbf{A}=\mathbf{V}\mathbf{D}\mathbf{V}^{-1}$, with eigenvalues $\lambda_1,\ldots,\lambda_n$. If $\mu$ is an eigenvalue of $\mathbf{A}+\mathbf{E}$ for a complex matrix $\mathbf{E}$, then

```{math}
:label: bauerfike
\min_{j=1,\ldots,n} |\mu - \lambda_j| \le \kappa(\mathbf{V}) \, \| \mathbf{E} \|\,,
```

where $\|\cdot\|$ and $\kappa$ are in the 2-norm.
````

The Bauer–Fike theorem tells us that eigenvalues can be perturbed by an amount that is $\kappa(\mathbf{V})$ times larger than perturbations to the matrix. This result is a bit less straightforward than it might seem—eigenvectors are not unique, so there are multiple possible values for $\kappa(\mathbf{V})$. Even so, the theorem indicates caution when a matrix has eigenvectors that form an ill-conditioned matrix. The limiting case of $\kappa(\mathbf{V})=\infty$ might be interpreted as indicating a nondiagonalizable matrix $\mathbf{A}$. The other extreme is also of interest: $\kappa(\mathbf{V})=1$, which implies that $\mathbf{V}$ is unitary.

```{index} ! normal matrix, unitary matrix
```

::::{prf:definition} Normal matrix
:label: definition-normalmatrix
If $\mathbf{A}$ has an EVD {eq}`evdecomp` with a unitary eigenvector matrix $\mathbf{V}$, then $\mathbf{A}$ is a {term}`normal matrix`.
::::

As we will see in {numref}`section-matrixanaly-symm-eig`, hermitian and real symmetric matrices are normal. Since the condition number of a unitary matrix is equal to 1, {eq}`bauerfike` guarantees that a perturbation of a normal matrix changes the eigenvalues by the same amount or less.

::::{prf:example} Eigenvalue conditioning
:label: demo-evd-bauerfike

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-evd-bauerfike-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-evd-bauerfike-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-evd-bauerfike-python
:::
```` 
`````

::::

## Computing the EVD

Roots of the characteristic polynomial are not used in numerical methods for finding eigenvalues.[^eigpoly] Practical algorithms for computing the EVD go beyond the scope of this book. The essence of the matter is the connection to matrix powers indicated in {eq}`evdpower`. (We will see much more about the importance of matrix powers in Chapter 8.)

If the eigenvalues have different complex magnitudes, then as $k\to\infty$ the entries on the diagonal of $\mathbf{D}^k$ become increasingly well separated and easy to pick out. It turns out that there is an astonishingly easy and elegant way to accomplish this separation without explicitly computing the matrix powers.

[^eigpoly]: In fact, the situation is reversed: eigenvalue methods are among the best ways to compute the roots of a given polynomial.

::::{prf:example} Francis QR iteration
:label: demo-evd-francisqr

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-evd-francisqr-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-evd-francisqr-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-evd-francisqr-python
:::
```` 
`````

::::

```{index} ! Francis QR iteration
```

The process demonstrated in @demo-evd-francisqr is known as the *Francis QR iteration*, and it can be formulated as an $O(n^3)$ algorithm for finding the EVD. It forms the basis of most practical eigenvalue computations, at least until the matrix size approaches $10^4$ or so.

## Exercises

``````{exercise}
:label: problem-evd-norm
**(a)** ✍ Suppose that matrix $\mathbf{A}$ has an eigenvalue $\lambda$. Show that for any induced matrix norm, $\| \mathbf{A} \|\ge |\lambda|$.

**(b)** ✍ Find a matrix $\mathbf{A}$ such that $\| \mathbf{A} \|_2$ is strictly larger than $|\lambda|$ for all eigenvalues $\lambda$. (Proof-by-computer isn't allowed here. You don't need to compute $\| \mathbf{A} \|_2$ exactly, just a lower bound for it.)
``````

``````{exercise}
:label: problem-evd-defective
✍ Prove that the matrix $\mathbf{B}$ in {eq}`jordanblock` does not have two independent eigenvectors.

``````
<!-- TODO: Get rid of rank. -->

``````{exercise}
:label: problem-evd-residual
⌨ In each part, find all the eigenvalues of $\mathbf{A}$. Then, choose one eigenvalue $\lambda$ and associated eigenvector $\mathbf{v}$ and compute $\twonorm{\mathbf{A} \mathbf{v} - \lambda \mathbf{v}}$, which should be comparable to machine epsilon.

**(a)** $\mathbf{A} = \begin{bmatrix}
2  & -1 & 0 \\
-1 &  2 & -1 \\
0  & -1 & 2
\end{bmatrix}$

**(b)** $\mathbf{A} = \begin{bmatrix}
2 & -1 & -1 \\
-2 &  2 & -1 \\
-1 & -2 & 2
\end{bmatrix}$

**(c)** $ \mathbf{A} = \begin{bmatrix}
2 & -1 & -1 \\
-1 &  2 & -1 \\
-1 & -1 & 2
\end{bmatrix} $

**(d)** $\mathbf{A} = \begin{bmatrix}
3 & 1 & 0 & 0 \\
1 & 3 & 1 & 0 \\
0 & 1 & 3 & 1 \\
0 & 0 & 1 & 3
\end{bmatrix}\qquad $

**(e)** $\mathbf{A} = \begin{bmatrix}
4 & -3 & -2 & -1\\
-2 &  4 & -2 & -1 \\
-1 & -2 & 4  & -1 \\
-1 & -2 & -1 & 4 \\
\end{bmatrix} $
``````

``````{exercise}
:label: problem-evd-triangular
 **(a)** ✍ Show that the eigenvalues of a diagonal $n\times n$ matrix $\mathbf{D}$ are the diagonal entries of $\mathbf{D}$. (That is, produce the associated eigenvectors.)

**(b)** ✍ The eigenvalues of a triangular matrix are its diagonal entries. Prove this in the $3\times 3$ case,

```{math}
\mathbf{T} =
\begin{bmatrix}
t_{11} & t_{12}&  t_{13}\\ 0 & t_{22} & t_{23} \\ 0 & 0 & t_{33}
\end{bmatrix},
```

by finding the eigenvectors. (Start by showing that $[1,0,0]^T$ is an eigenvector. Then show how to make $[a,1,0]^T$ an eigenvector, except for one case that does not change the outcome. Continue the same logic for $[a,b,1]^T$.)
``````

``````{exercise}
:label: problem-evd-matrixpoly
✍ Let $\mathbf{A}=\displaystyle\frac{\pi}{6}\begin{bmatrix} 4 & 1 \\ 4 & 4 \end{bmatrix}$.

**(a)** Show that

$$
\lambda_1=\pi,\, \mathbf{v}_1=\begin{bmatrix}1 \\ 2 \end{bmatrix}, \quad \lambda_2=\frac{\pi}{3},\, \mathbf{v}_2=\begin{bmatrix}1 \\ -2 \end{bmatrix}
$$

yield an EVD of $\mathbf{A}$.

**(b)** Use {eq}`matrixpolyevd` to evaluate $p(\mathbf{A})$, where $p(x) = (x-\pi)^4$.

**(c)** Use the function analog of {eq}`matrixpolyevd` to evaluate $\cos(\mathbf{A})$.
``````

``````{exercise}
:label: problem-evd-lumpstring
⌨ In @problem-linearsystems-lumpstring, you showed that the
displacements of point masses placed along a string satisfy a linear system $\mathbf{A}\mathbf{q}=\mathbf{f}$ for an $(n-1)\times(n-1)$ matrix $\mathbf{A}$. The eigenvalues and eigenvectors of $\mathbf{A}$ correspond to resonant frequencies and modes of vibration of the string. For $n=40$ and the physical parameters given in part (b) of that exercise, find the eigenvalue decomposition of $\mathbf{A}$. Report the three eigenvalues with smallest absolute value, and plot all three associated eigenvectors on a single graph (as functions of the vector row index).
``````

``````{exercise}
:label: problem-evd-francisqr
⌨ @demo-evd-francisqr suggests that the result of the Francis QR iteration as $k\to\infty$ sorts the eigenvalues on the diagonal according to a particular ordering. Following the code there as a model, create a random matrix with eigenvalues equal to $-9.6,-8.6,\ldots,10.4$, perform the iteration 200 times, and check whether the sorting criterion holds in your experiment as well.
``````

``````{exercise}
:label: problem-evd-random
⌨ Eigenvalues of random matrices and their perturbations can be very interesting.

**(a)** Let `A=randn(60,60)`.[^randn] Scatter plot its eigenvalues in the complex plane, using a plot aspect ratio of 1 and red diamonds as markers.

**(b)** Let $\mathbf{E}$ be another random $60\times 60$ matrix, and on top of the previous graph, plot the eigenvalues of $\mathbf{A}+0.05\mathbf{E}$ as blue dots. Repeat this for 100 different values of $\mathbf{E}$.

**(c)** Let `T=triu(A)`. On a new graph, scatter plot the eigenvalues of $\mathbf{T}$ in the complex plane. (They all lie on the real axis.)

**(d)** Repeat part (b) with $\mathbf{T}$ in place of $\mathbf{A}$.

**(e)** Compute some condition numbers and apply @theorem-bauer-fike to explain the dramatic difference between your plots with respect to the dot distributions.
``````

``````{exercise}
:label: problem-evd-functions
Suppose that $\mathbf{A}$ is diagonalizable and that @matrixpolyevd is used to define $\cos(\mathbf{A})$ and $\sin(\mathbf{A})$, with those functions substituted in for $p$ in the equation. Is it necessarily true that $\cos(\mathbf{A})^2+\sin(\mathbf{A})^2$ is an identity matrix? Explain why or why not.
``````

[^randn]: The `randn` function generates random numbers from a standard normal distribution. In Python, it is found in the `numpy.random` module.
