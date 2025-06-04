---
numbering:
  enumerator: 2.9.%s
---
(section-linsys-structure)=
# Exploiting matrix structure

A common situation in computation is that a problem has certain properties or structure that can be used to get a faster or more accurate solution. There are many properties of a matrix that can affect LU factorization. For example, an $n \times n$ matrix $A$ is **diagonally dominant** if

```{math}
:label: diag-dominant
  |A_{ii}| > \sum_{\substack{j=1\\ j \neq i}}^{n} |A_{ij}| \hskip 0.25in \text{for each } i=1,\ldots,n.
```

It turns out that a diagonally dominant matrix is guaranteed to be invertible, and row pivoting is not required for stability.

We next consider three important types of matrices that cause the LU factorization to be specialized in some important way.

## Banded matrices

```{index} ! bandwidth of a matrix, ! tridiagonal matrix
```

::::{prf:definition} Bandwidth
:label: definition-bandwidth
A matrix $\mathbf{A}$ has **upper bandwidth** $b_u$ if $j-i > b_u$ implies $A_{ij}=0$, and **lower bandwidth** $b_\ell$ if $i-j > b_\ell$ implies $A_{ij}=0$. We say the total {term}`bandwidth` is $b_u+b_\ell+1$. When $b_u=b_\ell=1$, we have the important case of a **tridiagonal matrix**. 
::::

(demo-structure-banded)=
::::{prf:example} Banded matrices
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #demo-structure-banded-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-structure-banded-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-structure-banded-python
:::
```` 
`````
::::


```{index} pivoting
```

If row pivoting is not used, the $\mathbf{L}$ and $\mathbf{U}$ factors preserve the lower and upper bandwidths of $\mathbf{A}$. This fact implies computational savings in both the factorization and the triangular substitutions because the zeros appear predictably and we can skip multiplication and addition with them. 

```{prf:observation}
The number of flops needed by LU factorization without pivoting is $O(b_u b_\ell n)$ when the upper and lower bandwidths are $b_u$ and $b_\ell$.
```

```{index} sparse matrix
```

In order to exploit the savings offered by sparsity, we would need to make modifications to {numref}`Function %s <function-lufact>` and the triangular substitution routines. Alternatively, we can get Julia to take advantage of the structure automatically by converting the matrix into a special type called **sparse**. Sparse matrices are covered in more detail in Chapters 7 and 8.

## Symmetric matrices

```{index} symmetric matrix
```

::::{prf:definition} Symmetric matrix
A square matrix $\mathbf{A}$ satisfying $\mathbf{A}^T = \mathbf{A}$ is called **symmetric**.
::::

Symmetric matrices arise frequently in applications because many types of interactions, such as gravitation and social-network befriending, are inherently symmetric. Symmetry in linear algebra simplifies many properties and algorithms. As a rule of thumb, if your matrix has symmetry, you want to exploit and preserve it. 

In $\mathbf{A}=\mathbf{L}\mathbf{U}$ we arbitrarily required the diagonal elements of $\mathbf{L}$, but not $\mathbf{U}$, to be one. That breaks symmetry, so we need to modify the goal to

$$
\mathbf{A}=\mathbf{L}\mathbf{D}\mathbf{L}^T,
$$

where $\mathbf{L}$ is unit lower triangular and $\mathbf{D}$ is diagonal. To find an algorithm for this factorization, we begin by generalizing {eq}`matrixouter` a bit without furnishing proof.

::::{prf:theorem} Linear combination of outer products
Let $\mathbf{D}$ be an $n\times n$ diagonal matrix with diagonal elements $d_1,d_2,\ldots,d_n$, and suppose $\mathbf{A}$ and $\mathbf{B}$ are $n\times n$ as well. Write the columns of $\mathbf{A}$ as $\mathbf{a}_1,\dots,\mathbf{a}_n$ and the rows of $\mathbf{B}$ as $\mathbf{b}_1^T,\dots,\mathbf{b}_n^T$. Then

```{math}
:label: matrixouter3
\mathbf{A}\mathbf{D}\mathbf{B} = \sum_{k=1}^n d_k \mathbf{a}_k \mathbf{b}_k^T.
```
::::

Let's derive the LDL$^T$ factorization for a small example.

(demo-structure-symm)=
::::{prf:example} Symmetric LDL$^T$ factorization
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #demo-structure-symm-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-structure-symm-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-structure-symm-python
:::
```` 
`````
::::

In practice we don't actually have to carry out any arithmetic in the upper triangle of $\mathbf{A}$ as we work, since the operations are always the mirror image of those in the lower triangle. As a result, it can be shown that LDL$^T$ factorization takes about half as much work as the standard LU.

::::{prf:observation}
LDL$^T$ factorization on an $n \times n$ symmetric matrix, when successful, takes $\sim \frac{1}{3}n^3$ flops.
::::

Just as pivoting is necessary to stabilize LU factorization, the LDL$^T$ factorization without pivoting may be unstable or even fail to exist. We won't go into the details, because our interest is in specializing the factorization to matrices that also possess another important property.

(sec-SPD)=
## Symmetric positive definite matrices

Suppose that $\mathbf{A}$ is $n\times n$ and $\mathbf{x}\in\mathbb{R}^n$. Observe that $\mathbf{x}^T\mathbf{A}\mathbf{x}$ is the product of $1\times n$, $n\times n$, and $n\times 1$ matrices, so it is a scalar, sometimes referred to as a **quadratic form**. It can be expressed as

```{math}
:label: quadratic-form
  \mathbf{x}^T\mathbf{A}\mathbf{x} = \sum_{i=1}^n \sum_{j=1}^n A_{ij}x_ix_j.
```

```{index} ! symmetric positive definite matrix
```

```{index} see: SPD matrix; symmetric positive definite matrix
```

::::{prf:definition} Symmetric positive definite matrix
A real $n\times n$ matrix $\mathbf{A}$ is called a **symmetric positive definite matrix** (or SPD matrix) if it is symmetric and, for all nonzero $\mathbf{x}\in\mathbb{R}^n$,

```{math}
:label: SPD-def
  \mathbf{x}^T \mathbf{A} \mathbf{x} > 0.
```
::::

The definiteness property is usually difficult to check directly from the definition. There are some equivalent conditions, though. For instance, a symmetric matrix is positive definite if and only if its eigenvalues are all real positive numbers. SPD matrices have important properties and appear in applications in which the definiteness is known for theoretical reasons.

Let us consider what definiteness means to the LDL$^T$ factorization. We compute

```{math}
  0 < \mathbf{x}^T\mathbf{A}\mathbf{x} = \mathbf{x}^T \mathbf{L} \mathbf{D} \mathbf{L}^T \mathbf{x} = \mathbf{z}^T \mathbf{D} \mathbf{z},
```

where $\mathbf{z}=\mathbf{L}^T \mathbf{x}$. Note that since $\mathbf{L}$ is unit lower triangular, it is nonsingular, so $\mathbf{x}=\mathbf{L}^{-T}\mathbf{z}$. By taking $\mathbf{z}=\mathbf{e}_k$ for $k=1,\ldots,n$, we can read the equalities from right to left and conclude that $D_{kk}>0$ for all $k$. That permits us to write a kind of square root formula:[^sqrt]

[^sqrt]: Except for this diagonal, positive definite case, it's not trivial to define the square root of a matrix, so don't generalize the notation used here.

```{math}
:label: diag-sqrt
  \mathbf{D} =
  \begin{bmatrix}
    D_{11} &        &        & \\
           & D_{22} &        & \\
           &        & \ddots & \\
           &        &        & D_{nn}
  \end{bmatrix}
=   \begin{bmatrix}
    \sqrt{D_{11}} &        &        & \\
           & \sqrt{D_{22}} &        & \\
           &        & \ddots & \\
           &        &        & \sqrt{D_{nn}}
  \end{bmatrix}^{\,2}
= \bigl( \mathbf{D}^{1/2} \bigr)^2.
```

Now we have $\mathbf{A}=\mathbf{L}\mathbf{D}^{1/2}\mathbf{D}^{1/2}\mathbf{L}^T= \mathbf{R}^T \mathbf{R}$, where $\mathbf{R} =\mathbf{D}^{1/2}\mathbf{L}^T$ is an upper triangular matrix whose diagonal entries are positive.

```{index} ! matrix factorization; Cholesky, Cholesky factorization
```

```{index} pivoting
```

::::{prf:theorem} Cholesky factorization
Any SPD matrix $\mathbf{A}$ may be factored as 

$$
\mathbf{A} = \mathbf{R}^T \mathbf{R},
$$

where $\mathbf{R}$ is an upper triangular matrix with positive diagonal elements. This is called the **Cholesky factorization**.
::::

While the unpivoted LDL$^T$ factorization is not stable and not even always possible, in the SPD case one can prove that pivoting is not necessary for the existence nor the stability of the Cholesky factorization. 

::::{prf:observation}
Cholesky factorization of an $n \times n$ SPD matrix takes $\sim \frac{1}{3}n^3$ flops.
::::

The speed and stability of the Cholesky factorization make it the top choice for solving linear systems with SPD matrices. As a side benefit, the Cholesky algorithm fails (in the form of an imaginary square root or division by zero) if and only if the matrix $\mathbf{A}$ is not positive definite. This is often the best way to test the definiteness of a symmetric matrix about which nothing else is known.

(demo-structure-cholesky)=
::::{prf:example} Cholesky factorization
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #demo-structure-cholesky-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-structure-cholesky-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-structure-cholesky-python
:::
```` 
`````
::::

## Exercises

``````{exercise}
:label: problem-structure-diagdom
✍  For each matrix, use {eq}`diag-dominant` to determine whether it is diagonally dominant.

```{math}
\mathbf{A} =
\begin{bmatrix}
3  & 1  & 0 & 1  \\
0  & -2 & 0 & 1  \\
-1 & 0  & 4 & -1 \\
0  & 0  & 0 & 6
\end{bmatrix},
\quad
\mathbf{B} =
\begin{bmatrix}
1  & 0  & 0  & 0 & 0  \\
0  & 1  & 0  & 0 & 0  \\
0  & 0  & 1  & 0 & 0  \\
0  & 0  & 0  & 1 & 0  \\
0  & 0  & 0  & 0 & 0
\end{bmatrix},
\quad \mathbf{C} =
\begin{bmatrix}
2  & -1 & 0  & 0      \\
-1 & 2  & -1 & 0      \\
0  & -1 & 2  & -1     \\
0  & 0  & -1 & 2
\end{bmatrix}.
```
``````

``````{exercise}
:label: problem-structure-SPD

⌨ For each matrix, use inspection or `cholesky` in Julia to determine whether it is SPD.

```{math}
\mathbf{A} =
\begin{bmatrix}
1 & 0 & -1 \\ 0 & 4 & 5 \\ -1 & 5 & 10
\end{bmatrix},
\qquad
\mathbf{B} =
\begin{bmatrix}
1 & 0 & 1 \\ 0 & 4 & 5 \\ -1 & 5 & 10
\end{bmatrix},
\qquad
\mathbf{C} =
\begin{bmatrix}
1 & 0 & 1 \\ 0 & 4 & 5 \\ 1 & 5 & 1
\end{bmatrix}.
```
``````

``````{exercise}
:label: problem-structure-SPDdiag

✍ Show that the diagonal entries of a symmetric positive definite matrix are positive numbers. (Hint: Apply certain special cases of {eq}`SPD-def`.)
``````

``````{exercise}
:label: problem-structure-luband

⌨ Using {numref}`Function {number} <function-lufact>` as a guide, write a function

``` julia
function luband(A,upper,lower)
```

that accepts upper and lower bandwidth values and returns LU factors (without pivoting) in a way that avoids doing arithmetic using the locations that are known to stay zero. (Hint: Refer to the more efficient form of `lufact` given in {numref}`section-linsys-efficiency`.)

Test your function on the matrix with elements

$$
A_{ij} = \begin{cases} \frac{1}{i+j}, & -1 \le i-j \le 2,\\ 
0 & \text{otherwise.} \end{cases}
$$
``````

``````{exercise}
:label: problem-structure-tridiagonal

⌨ The `Tridiagonal` matrix type invokes a specialized algorithm for solving a linear system. 

**(a)** Set `n=1000` and `t=0`.  In a loop that runs 50 times, generate a linear system via 
``` julia
A = triu( tril(rand(n,n),1), -1)
b = ones(n)
```
Using `@elapsed`, increment `t` by the time it takes to perform `A\b`. Print out the final value of `t`.

**(b)** Repeat the experiment of part (a), but generate the matrix via
``` julia
A = Tridiagonal(rand(n,n))
```
What is the ratio of running times for part (a) and (b)?

**(c)** Now perform the experiment of part (b) for $n=1000,1200,1400,\ldots,3000$, keeping the total time for each value of $n$ in a vector. Plot running time as a function of $n$ on a log-log scale. Is the time most like $O(n)$, $O(n^2)$, or $O(n^3)$? (If the answer is unclear, try increasing the number of solves per value of $n$ to 100 or more.)
``````

``````{exercise}
:label: problem-structure-ATA
✍ Prove that if $\mathbf{A}$ is any real invertible square matrix, then $\mathbf{A}^T\mathbf{A}$ is SPD. (Hint: First show that $\mathbf{x}^T\mathbf{A}^T\mathbf{A}\mathbf{x} \ge 0$ for all $\mathbf{x}$. Then explain why zero is ruled out if $\mathbf{x}\neq \boldsymbol{0}$.)

``````
