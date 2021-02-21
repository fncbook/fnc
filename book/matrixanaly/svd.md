# Singular value decomposition

We now introduce another factorization that is as fundamental as the EVD.

```{index} unitary matrix
```

````{proof:theorem}
  Let $\mathbf{A}\in\mathbb{C}^{m\times n}$. Then $\mathbf{A}$ can be written as

```{math}
    :label: svd
    \mathbf{A} = \mathbf{U} \mathbf{S} \mathbf{V}^*,
```

where $\mathbf{U}\in\mathbb{C}^{m\times m}$ and $\mathbf{V}\in\mathbb{C}^{n\times n}$ are unitary and $\mathbf{S}\in\mathbb{R}^{m\times n}$ is real and diagonal with nonnegative entries. If $\mathbf{A}$ is real, then so are $\mathbf{U}$ and $\mathbf{V}$ (which are then orthogonal matrices).
````

```{index} matrix; factorization
```
```{index} singular value decomposition
```
```{index} singular value
```

Equation {eq}`svd` is called a {term}`singular value decomposition`, or SVD, of $\mathbf{A}$. The columns of $\mathbf{U}$ and $\mathbf{V}$ are called **left** and **right singular vectors**, respectively. The diagonal entries of $\mathbf{S}$, written $\sigma_1,\ldots,\sigma_r$, for $r=\min\{m,n\}$, are called the **singular values** of $\mathbf{A}$. By convention the singular values are ordered so that

```{math}
:label: svdorder
\sigma_1 \ge \sigma_2 \ge \cdots \ge \sigma_r\ge 0, \qquad r=\min\{m,n\}.
```

We call $\sigma_1$ the **principal singular value** and $\mathbf{u}_{1}$ and $\mathbf{v}_{1}$ the **principal singular vectors**. The matrix $\mathbf{S}$ in the SVD is uniquely defined when the ordering is imposed, but the singular vectors are not---one could replace both $\mathbf{U}$ and $\mathbf{V}$ by their negatives, for example.


````{proof:example}
  Suppose $\mathbf{A}$ is a real matrix and that $\mathbf{A}=\mathbf{U}\mathbf{S}\mathbf{V}^T$ is an SVD. Then $\mathbf{A}^T=\mathbf{V}\mathbf{S}^T\mathbf{U}^T$ meets all the requirements of an SVD for $\mathbf{A}^T$: the first and last matrices are orthogonal, and the middle matrix is diagonal with nonnegative entries. Hence $\mathbf{A}$ and $\mathbf{A}^T$ have the same singular values.
````


## Interpreting the SVD

Another way to write $\mathbf{A}=\mathbf{U}\mathbf{S}\mathbf{V}^*$ is as $\mathbf{A}\mathbf{V}=\mathbf{U}\mathbf{S}$. Taken columnwise, this means

```{math}
:label: svdcolumns
\mathbf{A} \mathbf{v}_{k} = \sigma_k \mathbf{u}_{k}, \qquad k=1,\ldots,r=\min\{m,n\}.
```

In words, each right singular vector is mapped by $\mathbf{A}$ to a scaled version of its corresponding left singular vector; the magnitude of scaling is its singular value.

Both the SVD and the EVD describe a matrix in terms of some special vectors and a small number of scalars. {ref}`tab-evdsvd` summarizes the key differences. The SVD sacrifices having the same basis in both source and image spaces---after all, they may not even have the same dimension---but as a result gains orthogonality in both spaces.

:::{list-table} Comparison of the EVD and SVD
:header-rows: 1
:name: tab-evdsvd

* - EVD
  - SVD
* - exists for most square matrices
  - exists for all rectangular and square matrices 
* - $\mathbf{A}\mathbf{x}_k = \lambda_k \mathbf{x}_k$ 
  - $\mathbf{A} \mathbf{v}_k = \sigma_k \mathbf{u}_k$ 
* - same basis for domain and range of $\mathbf{A}$ 
  - two orthogonal bases 
* - may have poor conditioning 
  - perfectly conditioned 
:::

## SVD and the 2-norm

The SVD is intimately connected to the 2-norm,} as the following theorem describes. 
```{index} norm; matrix
```

::::{proof:theorem} SVD properties
Let $\mathbf{A}\in\mathbb{C}^{m\times n}$ have an SVD $\mathbf{A}=\mathbf{U}\mathbf{S}\mathbf{V}^*$ in
which {eq}`svdorder` holds. Then:

1. The 2-norm satisfies

```{math}
:label: svdnorm
\| \mathbf{A} \|_2 = \sigma_1.
```

2. The rank of $\mathbf{A}$ is the number of nonzero singular values.
\item Let $r=\min\{m,n\}$. Then

```{math}
:label: svdcond
\kappa_2(\mathbf{A}) = \|\mathbf{A}\|_2\|\mathbf{A}^+\|_2 = \frac{\sigma_1}{\sigma_r},
```

where a division by zero implies that $\mathbf{A}$ does not have full rank.
::::


The conclusion {eq}`svdnorm` is can be proved by straightforward vector calculus (see {ref}`prob-svd-svdnormproof`). In the square case $m=n$, $\mathbf{A}$ having full rank is identical to being nonsingular. The SVD is the usual means for computing the 2-norm and condition number of a matrix. 

:::{sidebar} Demo
:class: demo
{doc}`demos/svd-props`
:::

## Connections to the EVD

```{index} matrix; factorization
```
```{index} matrix; hermitian
```

Let $\mathbf{A}=\mathbf{U}\mathbf{S}\mathbf{V}^*$ be $m\times n$, and consider the square 
 hermitian matrix $\mathbf{B}=\mathbf{A}^*\mathbf{A}$:

```{math}
\mathbf{B} = (\mathbf{V}\mathbf{S}^*\mathbf{U}^*) (\mathbf{U}\mathbf{S}\mathbf{V}^*) = \mathbf{V}\mathbf{S}^*\mathbf{S}\mathbf{V}^* = \mathbf{V}(\mathbf{S}^T\mathbf{S})\mathbf{V}^{-1}.
```

Note that $\mathbf{S}^T\mathbf{S}$ is a diagonal $n \times n$ matrix. There are two cases to consider. If $m \ge n$, then 

$$
\mathbf{S}^T\mathbf{S} = 
\begin{bmatrix}
  \sigma_1^2 & & \\
  & \ddots & \\
  & & \sigma_n^2
\end{bmatrix}.
$$

On the other hand, if $m<n$, then 

$$
\mathbf{S}^T\mathbf{S} =
\begin{bmatrix}
  \sigma_1^2 & & & \\
  & \ddots & & \\
  & & \sigma_m^2 & \\
  & & & \boldsymbol{0}
\end{bmatrix}.
$$

(The lower-right zero in the last matrix is $n-m$ square.) In both cases we may conclude that the squares of the singular values of $\mathbf{A}$ are all eigenvalues of $\mathbf{B}$. Conversely, an EVD of $\mathbf{B}$ reveals the singular values and a set of right singular vectors of $\mathbf{A}$. The left singular vectors could then be deduced from the identity $\mathbf{A}\mathbf{V} = \mathbf{U}\mathbf{S}$.

Another close connection between EVD and SVD comes via the $(m+n)\times (m+n)$ matrix

```{math}
:label: svdaugment
\mathbf{C} =
\begin{bmatrix}
0 & \mathbf{A}^* \\ \mathbf{A} & 0
\end{bmatrix}.
```

If $\sigma$ is a singular value of $\mathbf{B}$, then $\sigma$ and $-\sigma$ are eigenvalues of $\mathbf{C}$, and the associated eigenvector immediately reveals a left and a right singular vector (see {ref}`prob-svd-svdtoevd`). This connection is implicitly exploited by software to compute the SVD.

## Thin form

```{index} singular value decomposition; thin form
```

In [an earlier section](../leastsq/qr.md) we saw that a matrix has both a "full" and a "thin" form of the QR factorization. A similar situation holds with the SVD. Suppose $\mathbf{A}$ is $m\times n$ with $m > n$ and let $\mathbf{A}=\mathbf{U}\mathbf{S}\mathbf{V}^*$ be an SVD. The last $m-n$ rows of $\mathbf{S}$ are all zero due to the fact that $\mathbf{S}$ is diagonal. Hence

\begin{align*}
  \mathbf{U} \mathbf{S} & =
  \begin{bmatrix}
    \mathbf{u}_1 & \cdots & \mathbf{u}_n & \mathbf{u}_{n+1} & \cdots & \mathbf{u}_m
  \end{bmatrix}
  \begin{bmatrix}
    \sigma_1 & &  \\
    & \ddots &  \\
    & & \sigma_n \\
    & & \\
    & \boldsymbol{0} & \\
    & &
  \end{bmatrix} \notag\\
  &=
  \begin{bmatrix}
    \mathbf{u}_1 & \cdots & \mathbf{u}_n
  \end{bmatrix}
  \begin{bmatrix}
    \sigma_1 & &  \\
    & \ddots &  \\
    & & \sigma_n
  \end{bmatrix} = \widehat{\mathbf{U}} \widehat{\mathbf{S}},
\end{align*}
<!-- :label: svd-reduced -->

```{index} ONC
```

in which $\widehat{\mathbf{U}}$ is $m\times n$ and $\widehat{\mathbf{S}}$ is $n\times n$. This allows us to define the **thin SVD** $\mathbf{A}=\widehat{\mathbf{U}}\widehat{\mathbf{S}}\mathbf{V}^*$, in which $\widehat{\mathbf{S}}$ is square and diagonal and $\widehat{\mathbf{U}}$ is ONC but not unitary. This form is computationally preferable when $m\gg n$, since it requires far less storage and contains the same information about $\mathbf{A}$. 

<!-- 
\begin{exercises}
	\input{matrixanaly/exercises/SVD}
\end{exercises} -->

