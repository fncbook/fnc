# Symmetry and definiteness

```{index} matrix; symmetric
```

As we saw [earlier](../linsys/structure.md), symmetry can simplify the LU factorization into a symmetric form, $\mathbf{A}=\mathbf{L}\mathbf{D}\mathbf{L}^T$. Certain specializations occur as well for the eigenvalue and singular value factorizations. In this section we stay with complex-valued matrices, so we are interested in the case when $\mathbf{A}^*=\mathbf{A}$, or $\mathbf{A}$ is hermitian. However, we often loosely speak of "symmetry" to mean this property in the complex case. All of the statements in this section easily specialize to the real case.

## Unitary diagonalization

Suppose now that $\mathbf{A}^*=\mathbf{A}$ and that $\mathbf{A}=\mathbf{U}\mathbf{S}\mathbf{V}^*$ is an SVD. Since $\mathbf{S}$ is real and square, we have

```{math}
\mathbf{A}^* = \mathbf{V} \mathbf{S}^* \mathbf{U}^* = \mathbf{V} \mathbf{S} \mathbf{U}^*,
```

and it's tempting to conclude that $\mathbf{U}=\mathbf{V}$. Happily, this is nearly true. The following theorem is typically proved in an advanced linear algebra course.

(thm-spec-decomp)=
````{prf:theorem} Spectral decomposition
If $\mathbf{A}=\mathbf{A}^*$, then $\mathbf{A}$ has a diagonalization $\mathbf{A}=\mathbf{V} \mathbf{D} \mathbf{V}^{-1}$ in which $\mathbf{V}$ is unitary and $\mathbf{D}$ is diagonal and real.
````

Another way to state the result of this theorem is that a hermitian matrix has a complete set of orthonormal eigenvectors—that is, a *unitary diagonalization*—and real eigenvalues. In this case, the EVD 

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

is an SVD, where $|\mathbf{D}|$ is the elementwise absolute value and $\mathbf{T}$ is diagonal with $|T_{ii}|=1$ for all $i$.
````

::::{prf:proof}
Let $T_{ii}=\operatorname{sign}(D_{ii})$ for all $i$. Then $\mathbf{T}^2=\mathbf{I}$ and $|\mathbf{D}|=\mathbf{T}\mathbf{D}$. The result follows.
::::

The converse of the [Spectral Decomposition Theorem](thm-spec-decomp) is also true: every matrix with a unitary diagonalization and real eigenvalues is hermitian. However, there are non-hermitian matrices that meet just the requirement of a unitary EVD; any such matrix is called {term}`normal`.

::::{prf:example} Julia demo
:class: demo
:label: demos-symm-eig-normal
{doc}`demos/symm-eig-normal`
::::

Now consider again [the Bauer--Fike Theorem](thm-bauer-fike), which says that the condition number of the eigenvalues is bounded above by $\kappa(\mathbf{V})$, for any eigenvector matrix $\mathbf{V}$. Because $\kappa=1$ for any unitary or orthogonal matrix, the [Spectral Decomposition Theorem](thm-spec-decomp)  then implies that the condition number of the eigenvalues of a hermitian or any normal matrix is one. That is, eigenvalues of a normal matrix can be changed by no more than the norm of the perturbation to the matrix.

:::{prf:example} Julia demo
:class: demo
:label: demos-symm-eig-perturb
{doc}`demos/symm-eig-perturb`
:::

## Rayleigh quotient

```{index} Rayleigh quotient
```

Recall that for a matrix $\mathbf{A}$ and compatible vector $\mathbf{x}$, the quadratic form $\mathbf{x}^* \mathbf{A} \mathbf{x}$ is a scalar. With a suitable normalization, it becomes the {term}`Rayleigh quotient`

```{math}
:label: rayleigh
R_{\mathbf{A}}(\mathbf{x}) = \frac{ \mathbf{x}^* \mathbf{A} \mathbf{x}}{\mathbf{x}^* \mathbf{x}}.
```

If $\mathbf{v}$ is an eigenvector such that $\mathbf{A} \mathbf{v}=\lambda \mathbf{v}$, then one easily calculates that $R_{\mathbf{A}}(\mathbf{v})=\lambda.$ That is, the Rayleigh quotient maps an eigenvector into its associated eigenvalue.

If $\mathbf{A}^*=\mathbf{A}$, then the Rayleigh quotient has another interesting property: $\nabla R_{\mathbf{A}}(\mathbf{v})=\boldsymbol{0}$ if $\mathbf{v}$ is an eigenvector. By a multidimensional Taylor series, then,

```{math}
:label: rq-series
R_{\mathbf{A}}(\mathbf{v}+\epsilon\mathbf{z}) = R_{\mathbf{A}}(\mathbf{v}) + 0 + O( \epsilon^2) =  \lambda + O( \epsilon^2),
```

as $\epsilon\to 0$. The conclusion is that a good estimate of an eigenvector becomes an even better estimate of an eigenvalue.

:::{prf:example} Julia demo
:class: demo
:label: demos-symm-eig-rayleigh
{doc}`demos/symm-eig-rayleigh`
:::

## Definite and indefinite matrices


```{index} matrix; positive definite
```
```{index} hermitian positive definite
```

In the real case, we called a symmetric matrix $\mathbf{A}$ *symmetric positive definite* (SPD) if $\mathbf{x}^T \mathbf{A}\mathbf{x} > 0 $ for all nonzero vectors $\mathbf{x}$. In the complex case the relevant property is  {term}`hermitian positive definite` (HPD), meaning that $\mathbf{A}^*=\mathbf{A}$ and $\mathbf{x}^* \mathbf{A}\mathbf{x} > 0$ for all complex vectors $\mathbf{x}$. Putting this property together with the Rayleigh quotient leads to the following.

(thm-hpd)=
````{prf:theorem}
If $\mathbf{A}^*=\mathbf{A}$, then the following statements are equivalent.

1. $\mathbf{A}$ is HPD.
2. The eigenvalues of $\mathbf{A}$ are positive numbers.
3. Any unitary EVD of $\mathbf{A}$ is also an SVD of $\mathbf{A}$.
````

A hermitian matrix with all negative eigenvalues is called *negative definite*, and one with eigenvalues of different signs is *indefinite*. Finally, if one or more eigenvalues is zero and the rest have one sign, it is *semidefinite*.

## Exercises

1. ✍ Each line below is an EVD for a hermitian matrix. State whether the matrix is definite, indefinite, or semidefinite. Then state whether the given factorization is also an SVD, and if it is not, modify it to find an SVD.

    **(a)** $\begin{bmatrix}
      0 & 0 \\ 0 & -1
    \end{bmatrix} =   \begin{bmatrix}
      0 & 1 \\ 1 & 0
    \end{bmatrix}  \begin{bmatrix}
      -1 & 0 \\ 0 & 0
    \end{bmatrix}  \begin{bmatrix}
      0 & 1 \\ 1 & 0
    \end{bmatrix}, \qquad$
    **(b)** $\begin{bmatrix}
      4 & -2 \\ -2 & 1
    \end{bmatrix} = \begin{bmatrix}
      1 & -0.5 \\ -0.5 & -1
    \end{bmatrix}  \begin{bmatrix}
      5 & 0 \\ 0 & 0
    \end{bmatrix}  \begin{bmatrix}
      0.8 & -0.4 \\ -0.4 & -0.8
    \end{bmatrix}$
    <br><br>

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

2. ✍ Prove true, or give a counterexample: If $\mathbf{A}$ and $\mathbf{B}$ are hermitian matrices of the same size, then
   
    $$
    R_{\mathbf{A}+\mathbf{B}}(\mathbf{x}) = R_{\mathbf{A}}(\mathbf{x})+R_{\mathbf{B}}(x).
    $$

    ```{index} field of values
    ```
3. ⌨ The range of the function $R_{\mathbf{A}}(\mathbf{x})$ is a subset of the complex plane known as the *field of values* of the matrix $\mathbf{A}$. Use 500 random vectors to plot points in the field of values of $\mathbf{A} = \displaystyle  \begin{bmatrix}
     1  &   0   & -2\\
     0  &   2  &   0\\
    -2   &  0 &    1
    \end{bmatrix}$. Then compute its eigenvalues and guess what the exact field of values is.
    ::::{only} solutions
    It's the interval $[-1,3]$ (between the extreme eigenvalues).
    ::::

4. ✍ Let $\mathbf{A}=\displaystyle \begin{bmatrix}
    3 & -2 \\ -2 & 0
  \end{bmatrix}.$
  
    **(a)** Write out $R_{\mathbf{A}}(\mathbf{x})$ explicitly as a function of $x_1$ and $x_2$.

    **(b)** Find $R_{\mathbf{A}}(\mathbf{x})$ for $x_1=1$, $x_2=2$.
    
    **(c)** Find the gradient vector $\nabla R_{\mathbf{A}}(\mathbf{x})$.
    
    **(d)** Show that the gradient vector is zero when $x_1=1$, $x_2=2$.

5. ✍ A *skew-Hermitian* matrix is one that satisfies $\mathbf{A}^*=-\mathbf{A}$. Show that if $\mathbf{A}$ is skew-Hermitian, then $R_{\mathbf{A}}$ is imaginary-valued.

