---
numbering:
  enumerator: 7.4.%s
---
(section-matrixanaly-symm-eig)=

# Symmetry and definiteness

```{index} symmetric matrix
```

As we saw in {numref}`section-linsys-structure`, symmetry can simplify the LU factorization into the symmetric form $\mathbf{A}=\mathbf{L}\mathbf{D}\mathbf{L}^T$. Important specializations occur as well for the eigenvalue and singular value factorizations. In this section we stay with complex-valued matrices, so we are interested in the case when $\mathbf{A}^*=\mathbf{A}$, i.e., $\mathbf{A}$ is hermitian. However, we often loosely speak of *symmetry* to mean this property even in the complex case. All the statements in this section easily specialize to the real case.

## Normality

Suppose now that $\mathbf{A}^*=\mathbf{A}$ and that $\mathbf{A}=\mathbf{U}\mathbf{S}\mathbf{V}^*$ is an SVD. Since $\mathbf{S}$ is real and square, we have

```{math}
\mathbf{A}^* = \mathbf{V} \mathbf{S}^* \mathbf{U}^* = \mathbf{V} \mathbf{S} \mathbf{U}^*,
```

and it's tempting to conclude that $\mathbf{U}=\mathbf{V}$. Happily, this is nearly true. The following theorem is typically proved in an advanced linear algebra course.

````{prf:theorem} Spectral decomposition
:label: theorem-symm-eig-spectral
If $\mathbf{A}=\mathbf{A}^*$, then $\mathbf{A}$ has a diagonalization $\mathbf{A}=\mathbf{V} \mathbf{D} \mathbf{V}^{-1}$ in which $\mathbf{V}$ is unitary and $\mathbf{D}$ is diagonal and real.
````

Another way to state the result of this theorem is that a hermitian matrix has real eigenvalues and a complete set of orthonormal eigenvectors—that is, the matrix is normal. Because hermitian matrices are normal, their eigenvalue condition number is guaranteed to be 1 by {numref}`Theorem {number} <theorem-bauer-fike>`.

:::{note}
The converse of {numref}`Theorem {number} <theorem-symm-eig-spectral>` is also true: every normal matrix with real eigenvalues is hermitian. This was illustrated in @demo-evd-bauerfike.
:::

For a hermitian matrix, the EVD

$$
\mathbf{A}=\mathbf{V}\mathbf{D}\mathbf{V}^{-1}=\mathbf{V} \mathbf{D} \mathbf{V}^*
$$

is almost an SVD.

```{index} unitary matrix
```

````{prf:theorem}
If $\mathbf{A}^*=\mathbf{A}$ and $\mathbf{A}=\mathbf{V}\mathbf{D}\mathbf{V}^{-1}$ is a unitary diagonalization, then
  
```{math}
:label: herm-svd
\mathbf{A} = (\mathbf{V}\mathbf{T})\cdot |\mathbf{D}|\cdot \mathbf{V}^*
```

is an SVD, where $|\mathbf{D}|$ is the elementwise absolute value and $\mathbf{T}$ is diagonal with $|T_{ii}|=1$ for all $i$. In particular, the absolute values of the eigenvalues of $\mathbf{A}$ are the singular values of $\mathbf{A}$.
````

::::{prf:proof}
:enumerated: false

Let $T_{ii}=\operatorname{sign}(D_{ii})$ for all $i$. Then $\mathbf{T}^2=\mathbf{I}$, $|\mathbf{D}|=\mathbf{T}\mathbf{D}$, and

$$
\mathbf{A}=\mathbf{V} \mathbf{D} \mathbf{V}^*=\mathbf{V} \mathbf{T}^2 \mathbf{D} \mathbf{V}^*=(\mathbf{V} \mathbf{T}) (\mathbf{T} \mathbf{D}) \mathbf{V}^*.
$$
::::

## Rayleigh quotient

```{index} ! Rayleigh quotient
```

Recall that for a matrix $\mathbf{A}$ and compatible vector $\mathbf{x}$, the quadratic form $\mathbf{x}^* \mathbf{A} \mathbf{x}$ is a scalar.

::::{prf:definition}
:label: definition-rayleighquotient
Given hermitian $\mathbf{A}$ and nonzero vector $\mathbf{x}$, the {term}`Rayleigh quotient` is the function

```{math}
:label: rayleigh
R_{\mathbf{A}}(\mathbf{x}) = \frac{ \mathbf{x}^* \mathbf{A} \mathbf{x}}{\mathbf{x}^* \mathbf{x}}.
```

::::

If $\mathbf{v}$ is an eigenvector such that $\mathbf{A} \mathbf{v}=\lambda \mathbf{v}$, then one easily calculates that $R_{\mathbf{A}}(\mathbf{v})=\lambda.$ That is, the Rayleigh quotient maps an eigenvector into its associated eigenvalue.

If $\mathbf{A}^*=\mathbf{A}$, then the Rayleigh quotient has another interesting property: $\nabla R_{\mathbf{A}}(\mathbf{v})=\boldsymbol{0}$ if $\mathbf{v}$ is an eigenvector. By a multidimensional Taylor series, then,

```{math}
:label: rq-series
R_{\mathbf{A}}(\mathbf{v}+\epsilon\mathbf{z}) = R_{\mathbf{A}}(\mathbf{v}) + 0 + O( \epsilon^2) =  \lambda + O( \epsilon^2),
```

as $\epsilon\to 0$. The conclusion is that a good estimate of an eigenvector becomes an even better estimate of an eigenvalue.

::::{prf:example} Rayleigh quotient
:label: demo-symm-eig-rayleigh

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-symm-eig-rayleigh-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-symm-eig-rayleigh-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-symm-eig-rayleigh-python
:::
```` 
`````

::::

## Definite, semidefinite, and indefinite matrices

```{index} symmetric positive definite matrix
```

```{index} see: hermitian positive definite matrix; symmetric positive definite matrix
```

In the real case, we called a symmetric matrix $\mathbf{A}$ SPD if $\mathbf{x}^T \mathbf{A}\mathbf{x} > 0 $ for all nonzero vectors $\mathbf{x}$. There is an analogous definition for complex matrices.

::::{prf:definition} HPD matrix
:label: definition-HPDmatrix
An {term}`HPD matrix` is a complex $n\times n$ matrix $\mathbf{A}$ that is hermitian and for which

```{math}
:label: HPD-def
  \mathbf{x}^* \mathbf{A} \mathbf{x} > 0
```

for all nonzero $\mathbf{x}\in\mathbb{C}^n$.

::::

:::{warning}
It's pretty common to use the term *positive definite* without either the "symmetric" or "hermitian" adjective. This is confusing, because the definiteness property is of little value without the symmetry property, which is generally considered implied.
:::

Putting the HPD property together with the Rayleigh quotient leads to the following.

````{prf:theorem}
:label: theorem-symm-eig-hpd
If $\mathbf{A}^*=\mathbf{A}$, then the following statements are equivalent.

1. $\mathbf{A}$ is HPD.
2. The eigenvalues of $\mathbf{A}$ are positive numbers.
````

::::{prf:proof}
:enumerated: false

Suppose item 1 is true. If $\mathbf{A}\mathbf{x} = \lambda \mathbf{x}$ is an eigenpair, then a Rayleigh quotient implies that

$$
\lambda = \frac{ \mathbf{x}^*\mathbf{A}\mathbf{x} }{\mathbf{x}^*\mathbf{x}} > 0.
$$

Hence, item 2 is true. Conversely, suppose item 2 is known. Then we can write the EVD as $\mathbf{A}=\mathbf{V}\mathbf{S}^2\mathbf{V}^*$, where the $S_{ii}$ are positive square roots of the eigenvalues. Hence

$$
\mathbf{x}^*\mathbf{A}\mathbf{x} = \mathbf{x}^*\mathbf{V}\mathbf{S}^2\mathbf{V}^*\mathbf{x} = \|\mathbf{S}\mathbf{V}^*\mathbf{x}\|^2 > 0,
$$

as both $\mathbf{S}$ and $\mathbf{V}$ are invertible. Thus, item 1 is true.
::::

According to {numref}`Theorem {number} <theorem-symm-eig-hpd>`, for an HPD matrix, the EVD $\mathbf{A}=\mathbf{V}\mathbf{D}\mathbf{V}^*$ meets all the requirements of the SVD, provided the ordering of eigenvalues is chosen appropriately.

A hermitian matrix with all negative eigenvalues is called **negative definite**, and one with eigenvalues of different signs is **indefinite**. Finally, if one or more eigenvalues is zero and the rest have one sign, it is positive or negative **semidefinite**.

## Exercises

``````{exercise}
:label: problem-symmeig-definite
✍ Each line below is an EVD for a hermitian matrix. State whether the matrix is definite, indefinite, or semidefinite. Then state whether the given factorization is also an SVD, and if it is not, modify it to find an SVD.

**(a)** 
$\begin{bmatrix}
0 & 0 \\ 0 & -1
\end{bmatrix} =   \begin{bmatrix}
0 & 1 \\ 1 & 0
\end{bmatrix}  \begin{bmatrix}
-1 & 0 \\ 0 & 0
\end{bmatrix}  \begin{bmatrix}
0 & 1 \\ 1 & 0
\end{bmatrix}$

**(b)** 
$\begin{bmatrix}
4 & -2 \\ -2 & 1
\end{bmatrix} = \begin{bmatrix}
1 & -0.5 \\ -0.5 & -1
\end{bmatrix}  \begin{bmatrix}
5 & 0 \\ 0 & 0
\end{bmatrix}  \begin{bmatrix}
0.8 & -0.4 \\ -0.4 & -0.8
\end{bmatrix}$

**(c)**
$\begin{bmatrix}
-5 & 3\\ 3 & -5
\end{bmatrix} =  \begin{bmatrix}
\alpha & \alpha \\ \alpha & -\alpha
\end{bmatrix}  \begin{bmatrix}
-2 & 0 \\ 0 & -8
\end{bmatrix}  \begin{bmatrix}
\alpha & \alpha \\ \alpha & -\alpha
\end{bmatrix}, \quad\alpha=1/\sqrt{2}$
``````

``````{exercise}
:label: problem-symmeig-named
⌨ The matrix names below are found in `MatrixDepot` for Julia, `gallery` for MATLAB, and `rogues` for Python. You will have to adjust the syntax accordingly. For each matrix, determine whether it is positive definite, negative definite, positive or negative semidefinite, or indefinite. 

**(a)** `pei(5)` $ - 6 \mathbf{I}$

**(b)** `hilb(8)` $ - 2 \mathbf{I}$

**(c)** `dingdong(20)`

**(d)** `lehmer(100)`

**(e)** `fiedler(200)`
``````

``````{exercise}
:label: problem-symmeig-linearity
✍ Prove true, or give a counterexample: If $\mathbf{A}$ and $\mathbf{B}$ are hermitian matrices of the same size, then

$$
R_{\mathbf{A}+\mathbf{B}}(\mathbf{x}) = R_{\mathbf{A}}(\mathbf{x})+R_{\mathbf{B}}(\mathbf{x}).
$$
``````

```{index} field of values
```

``````{exercise}
:label: problem-symmeig-fov
⌨ The range of the function $R_{\mathbf{A}}(\mathbf{x})$ is a subset of the complex plane known as the *field of values* of the matrix $\mathbf{A}$. Use 500 random vectors to plot points in the field of values of $\mathbf{A} = \displaystyle  \begin{bmatrix}
1  &   0   & -2\\
0  &   2  &   0\\
-2   &  0 &    1
\end{bmatrix}$. Then compute its eigenvalues and guess what the exact field of values is.
``````

``````{exercise}
:label: problem-symmeig-gradient
✍ Let $\mathbf{A}=\displaystyle \begin{bmatrix} 3 & -2 \\ -2 & 0 \end{bmatrix}.$

**(a)** Write out $R_{\mathbf{A}}(\mathbf{x})$ explicitly as a function of $x_1$ and $x_2$.

**(b)** Find $R_{\mathbf{A}}(\mathbf{x})$ for $x_1=1$, $x_2=2$.

**(c)** Find the gradient vector $\nabla R_{\mathbf{A}}(\mathbf{x})$.

**(d)** Show that the gradient vector is zero when $x_1=1$, $x_2=2$.
``````

``````{exercise}
:label: problem-symmeig-skew
✍ A *skew-Hermitian* matrix is one that satisfies $\mathbf{A}^*=-\mathbf{A}$. Show that if $\mathbf{A}$ is skew-Hermitian, then $R_{\mathbf{A}}$ is imaginary-valued.
``````

``````{exercise}
:label: problem-symmeig-evd
⌨ Thanks largely to {numref}`Theorem {number} <theorem-symm-eig-spectral>`, the eigenvalue problem for symmetric/hermitian matrices is easier than for general matrices. 

**(a)** Let $\mathbf{A}$ be a $1000\times 1000$ random real matrix, and let $\mathbf{S}=\mathbf{A}+\mathbf{A}^T$. Time finding the eigenvalues of $\mathbf{A}$ and then of $\mathbf{S}$. You should find that the computation for $\mathbf{S}$ is around an order of magnitude faster.

**(b)** Perform the experiment from part (a) on $n\times n$ matrices for $n=200,300,\ldots,1600$. Plot running time as a function of $n$ for both matrices on a single log-log plot. Is the ratio of running times roughly constant, or does it grow with $n$?
``````
