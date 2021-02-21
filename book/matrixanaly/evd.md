# Eigenvalue decomposition

To this point we have dealt frequently with the solution of the linear system $\mathbf{A} \mathbf{A}\mathbf{x}=\mathbf{A}\mathbf{b}$. Alongside this problem in its importance to linear algebra is the eigenvalue problem,

```{math}
  :label: eigdef
  \mathbf{A} \mathbf{A}\mathbf{x} = \lambda \mathbf{A}\mathbf{x},
```

```{index} eigenvalue
```
```{index} eigenvector
```
for a scalar **eigenvalue** $\lambda$ and an associated nonzero **eigenvector** $\mathbf{x}$.

## Complex matrices

A matrix with real entries can have complex eigenvalues. Therefore we assume all matrices, vectors, and scalars may be complex in what follows. Recall that a complex number can be represented as $a+i b$, for real $a$ and $b$ and where $i^2=-1$. The **complex conjugate** of $x=a+i b$ is denoted $\bar{x}$ and is given by $\bar{x}=a-i b$. The {term}`hermitian` or conjugate transpose of a matrix $\mathbf{A}$ is denoted $\mathbf{A}^*$ and is given by $\mathbf{A}^*=(\overline{\mathbf{A}})^T=\overline{\mathbf{A}^T}$.

For the most part, "hermitian" replaces "transpose" when dealing with complex matrices. For instance, the inner product of complex vectors $\mathbf{A}\mathbf{u}$ and $\mathbf{A}\mathbf{v}$ is

```{math}
  :label: complexinnerprod
  \mathbf{A}\mathbf{u}^* \mathbf{A}\mathbf{v} = \sum_{k=1}^n \overline{u}_k v_k,
```

```{index} unitary matrix
```
```{index} matrix; unitary
```
```{index} orthogonal; matrix
```

which in turn defines the 2-norm for complex vectors, and thereby matrices as well. The definitions of orthogonal and orthonormal sets of complex-valued vectors use ${}^*$ instead of ${}^T$. The analog of an orthogonal matrix in the complex case---that is, a square matrix whose columns are orthonormal in the complex sense---is said to be {term}`unitary`. A unitary matrix $\mathbf{U}$ satisfies $\mathbf{U}^{-1}=\mathbf{U}^*$ and $\| \mathbf{U}\mathbf{x} \|_2=\| \mathbf{x} \|_2$ for any complex vector $\mathbf{x}\in\mathbb{C}^n$.


## Eigenvalue decomposition

```{index} characteristic polynomial
```

The eigenvalue equation $\mathbf{A}\mathbf{x}=\lambda\mathbf{x}$ is equivalent to $(\lambda\mathbf{I} - \mathbf{A})\mathbf{x}=\boldsymbol{0}$, which, since $\mathbf{x}$ is nonzero in order to be an eigenvector, implies that $\lambda\mathbf{I} - \mathbf{A}$ is a singular matrix. This observation leads to the familiar property of an eigenvalue being a root of the **characteristic polynomial** $\det(\lambda \mathbf{I} - \mathbf{A})$. From here one concludes that an $n\times n$ matrix has $n$ eigenvalues, counting multiplicity.

Hence suppose that $\mathbf{A}\mathbf{v}_k=\lambda_k\mathbf{v}_k$ for $k=1\ldots,n$. We can summarize these as

```{math}
:label: ev-all
\begin{align*}
   \begin{bmatrix}
    \mathbf{A}\mathbf{v}_1 & \mathbf{A}\mathbf{v}_2 & \cdots & \mathbf{A}\mathbf{v}_n
  \end{bmatrix}
  &=
    \begin{bmatrix}
      \lambda_1 \mathbf{v}_1 & \lambda_2\mathbf{v}_2 & \cdots & \lambda_n \mathbf{v}_n
    \end{bmatrix}\notag\\
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
  \end{bmatrix} \notag\\
  \mathbf{A} \mathbf{V} &= \mathbf{V} \mathbf{D}.
\end{align*}
```

So far $\mathbf{A}$ could be any square matrix. But if we also assume that $\mathbf{V}$ is a nonsingular matrix, we can rewrite {eq}`ev-all` as

```{math}
  :label: evdecomp
  \mathbf{A} = \mathbf{V} \mathbf{D} \mathbf{V}^{-1}.
```

```{index} matrix; factorization
```
```{index} eigenvalue decomposition
```
```{index} matrix; diagonalizable
```

Equation {eq}`evdecomp` is called an {term}`eigenvalue decomposition` (EVD) of $\mathbf{A}$. If $\mathbf{A}$ has an EVD, we say that $\mathbf{A}$ is {term}`diagonalizable`; otherwise $\mathbf{A}$ is **nondiagonalizable** (or *defective*). One simple example of a nondiagonalizable matrix is

```{math}
  :label: jordanblock
  \mathbf{B} = \begin{bmatrix}
    1 & 1\\0 & 1
  \end{bmatrix}.
```

There is a common circumstance in which we can guarantee an EVD exists; the proof of the following theorem can be found in many elementary texts on linear algebra.

````{proof:theorem}
If the $n\times n$ matrix $\mathbf{A}$ has $n$ distinct eigenvalues, then $\mathbf{A}$ is diagonalizable.
````

The {term}`eig` command can be used to compute the eigenvalue decomposition of a given matrix.

```{sidebar} Demo
:class: demo
{doc}`demos/evd-eigen`
```

observe that if $\mathbf{A}\mathbf{v} = \lambda \mathbf{v}$ for nonzero $\mathbf{v}$, then the equation remains true for any nonzero multiple of $\mathbf{v}$: eigenvectors are not unique, and neither is an EVD.

## Similarity and change of basis

```{index} matrix; similar
```

The particular relationship between matrices $\mathbf{A}$ and $\mathbf{D}$ in {eq}`evdecomp` is important. If $\mathbf{S}$ is any nonsingular matrix, we say that $\mathbf{B}=\mathbf{S}\mathbf{A}\mathbf{S}^{-1}$ is **similar** to $\mathbf{A}$. A similarity transformation does not change eigenvalues, a fact that is typically proved in elementary linear algebra texts.

````{proof:theorem}
If $\mathbf{X}$ is an nonsingular matrix, then $\mathbf{X}\mathbf{A}\mathbf{X}^{-1}$ has the same eigenvalues as $\mathbf{A}$.
````

Similarity transformation has a relatively simple interpretation. First, consider the product of a nonsingular $\mathbf{X}$ with any vector:

```{math}
  \mathbf{y} = \mathbf{X} \mathbf{z} = z_1 \mathbf{x}_1 +  \dots + z_n \mathbf{x}_n.
```

We call $z_1,\ldots,z_n$ the *coordinates* of the vector $\mathbf{y}$ with respect to the columns of $\mathbf{X}$. That is, $\mathbf{z}$ is a representation of $\mathbf{y}$ relative to the basis implied by the columns of $\mathbf{X}$. But then $\mathbf{z} = \mathbf{X}^{-1} \mathbf{y}$, so left-multiplication by $\mathbf{X}^{-1}$ converts the vector $\mathbf{y}$ into that representation. In other words, multiplication by the inverse of a matrix performs a **change of basis** into the coordinates associated with the matrix.

Now consider the EVD {eq}`evdecomp` and the product $\mathbf{u} =\mathbf{A} \mathbf{x}$, or $(\mathbf{V}^{-1}\mathbf{u}) = \mathbf{D}(\mathbf{V}^{-1}\mathbf{x})$. This equation says that if you express the input $\mathbf{x}$ and the output $\mathbf{u}$ into the coordinates of the $\mathbf{V}$-basis, then the relationship between them is diagonal. That is, the EVD is about finding a basis for $\mathbb{C}^n$ in which the map $\mathbf{x}\mapsto\mathbf{A}\mathbf{x}$ is a diagonal one in which the coordinates are independently rescaled.

The fact that the EVD represents a change of basis in both the domain and range spaces makes it useful for matrix powers: $\mathbf{A}^2=(\mathbf{V}\mathbf{D}\mathbf{V}^{-1})(\mathbf{V}\mathbf{D}\mathbf{V}^{-1})=\mathbf{V}\mathbf{D}(\mathbf{V}^{-1}\mathbf{V})\mathbf{D}\mathbf{V}^{-1}=\mathbf{V}\mathbf{D}^2\mathbf{V}^{-1}$, and so on. Because $\mathbf{D}$ is diagonal, its power $\mathbf{D}^k$ is just the diagonal matrix of the powers of the eigenvalues.

## Conditioning of eigenvalues


```{index} eigenvalue; conditioning of
```
```{index} condition number; of eigenvalues
```
```{index} Bauer--Fike theorem
```
[^BFT]: We will apply it only in the 2-norm, though it is more generally true.

Just as linear systems have condition numbers that quantify the effect of fixed precision, eigenvalue problems may be poorly conditioned too. While many possible results can be derived, we will use just one, the **Bauer--Fike theorem**.[^BFT]   

(thm-bauer-fike)=
````{proof:theorem} Bauer--Fike
Let $\mathbf{A}\in\mathbb{C}^{n\times n}$ be diagonalizable, $\mathbf{A}=\mathbf{V}\mathbf{D}\mathbf{V}^{-1}$, with eigenvalues $\lambda_1,\ldots,\lambda_n$. If $\mu$ is an eigenvalue of $\mathbf{A}+\mathbf{E}$ for a complex matrix $\mathbf{E}$, then

```{math}
:label: bauerfike
\min_{j=1,\ldots,n} |\mu - \lambda_j| \le \kappa(\mathbf{V}) \| \mathbf{E} \|,
```

where $\|\cdot\|$ and $\kappa$ are in the 2-norm.
````

The Bauer--Fike theorem tells us that eigenvalues can be perturbed by an amount that is $\kappa(\mathbf{V})$ times larger than perturbations to the matrix. This result is a bit less straightforward than it might seem---eigenvectors are not unique, so there are multiple possible values for $\kappa(\mathbf{V})$. Still, the theorem indicates caution when a matrix has eigenvectors that form an ill-conditioned matrix. The limiting case of $\kappa(\mathbf{V})=\infty$ might be interpreted as indicating a nondiagonalizable matrix $\mathbf{A}$.

```{index} matrix; normal
```

At the other extreme, if a unitary eigenvector matrix $\mathbf{V}$ can be found, then $\kappa(\mathbf{V})=1$ and {eq}`bauerfike` guarantees that eigenvalues are robust under perturbations to the original matrix $\mathbf{A}$. Such matrices are called **normal**, and they include the hermitian (or real symmetric) matrices. We consider them again [later on](symm-eig).

```{sidebar} Demo
:class: demo
{doc}`demos/evd-bauerfike`
```

## Computing the EVD

In elementary linear algebra you use the characteristic polynomial to compute the eigenvalues of small matrices. However, computing polynomial roots in finite time is impossible for degree 5 and over. In principle one could use Newton-like methods to find all of the roots, but doing so is relatively slow and difficult. Furthermore we know that polynomial roots tend to become poorly conditioned when roots get close to one another (see {ref}`example-quadrootcond`).

Practical algorithms for computing the EVD go beyond the scope of this book. The essence of the matter is the connection to matrix powers, that is, $\mathbf{A}^k = \mathbf{V} \mathbf{D}^k \mathbf{V}^{-1}$. (We will see much more about the importance of matrix powers in a [later chapter](../krylov).) If the eigenvalues have different complex magnitudes, then as $k\to\infty$ the entries on the diagonal of $\mathbf{D}^k$ become increasingly well separated and easy to pick out. It turns out that there is an astonishingly easy and elegant way to accomplish this separation without explicitly computing the matrix powers.

```{sidebar} Demo
:class: demo
{doc}`demos/evd-francisqr`
```

```{index} Francis QR iteration
```

The process demonstrated in {ref}`example-qriter` is known as the *Francis QR iteration*, and it can be formulated as an $O(n^3)$ algorithm for finding the EVD. Such an algorithm is what the `eigen` command uses.

<!-- 
\begin{exercises}
  \input{matrixanaly/exercises/EigenvalueDecomposition}
\end{exercises} -->
