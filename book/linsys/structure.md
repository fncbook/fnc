# Exploiting matrix structure

A common situation in computation is that a problem has certain properties or structure that can be used to get a faster or more accurate solution. There are many properties of a matrix that can affect LU factorization. For example, an $n \times n$ matrix $A$ is **diagonally dominant** if

```{math}
  :label: diag-dominant
  |A_{ii}| > \sum_{\substack{j=1\\ j \neq i}}^{n} |A_{ij}|, \hskip 0.25in \text{for each } i=1,\ldots,n.
```

It turns out that a diagonally dominant matrix is guaranteed to be nonsingular and row pivoting is not required for stability; i.e., $\mathbf{A}=\mathbf{L}\mathbf{U}$ is just as good as $\mathbf{P}\mathbf{A}=\mathbf{L}\mathbf{U}$.

We next consider three important types of matrices that cause the LU factorization to be specialized in some important way.

## Banded matrices

```{index} bandwidth
```

```{index} tridiagonal matrix
```

```{proof:example} Julia demo
:class: demo
{doc}`demos/structure-banded`
```

A matrix $\mathbf{A}$ has **upper bandwidth** $b_u$ if $j-i > b_u$ implies $A_{ij}=0$, and **lower bandwidth** $b_\ell$ if $i-j > b_\ell$ implies $A_{ij}=0$. We say the total {term}`bandwidth` is $b_u+b_\ell+1$. When $b_u=b_\ell=1$, we have the important case of a {term}`tridiagonal matrix`.  The {term}`spy` function in `Plots` is handy for visualizing the location of nonzero elements.

```{margin}
The number of flops needed by LU factorization is $O(b_ub_\ell n)$ when the upper and lower bandwidths are $b_u$ and $b_\ell$.
```

```{index} pivoting
```

If no row pivoting is used, the $\mathbf{L}$ and $\mathbf{U}$ factors preserve the lower and upper bandwidths of $\mathbf{A}$. This is a consequence of the row elimination process—fewer zeros need to be introduced into each column (equivalently, most of the row multipliers are zero), and adding rows downward cannot introduce new nonzeros into the upper triangle. This observation implies computational savings in both the factorization and the triangular substitutions, because the zeros appear predictably and we can skip multiplication and addition with them. To be precise, the number of flops needed by LU factorization is $O(b_ub_\ell n)$ when the upper and lower bandwidths are $b_u$ and $b_\ell$.

```{index} sparse matrix
```

```{proof:example} Julia demo
:class: demo
{doc}`demos/structure-timing`
```

In order to take advantage of the savings, we would need minor modifications to {ref}`function-lufact` and the triangular substitution routines. Alternatively, we can get Julia to take advantage of the structure automatically by converting the matrix into a special type called {term}`sparse`, defined in the `SparseArrays` package. Sparse matrices are covered in more detail in later chapters.

## Symmetric matrices

```{index} symmetric matrix
```

If $\mathbf{A}^T=\mathbf{A}$, then $\mathbf{A}$ is symmetric. Symmetric matrices arise frequently because so many types of interactions are symmetric: gravitation, social-network befriending, etc. Symmetry in linear algebra simplifies many properties and algorithms. As a rule of thumb, if your matrix has symmetry, you want to exploit and preserve it.

```{proof:example} Julia demo
:class: demo
{doc}`demos/structure-symm`
```

Now, if we create an LU factorization $\mathbf{A}=\mathbf{L}\mathbf{U}$ of a symmetric $\mathbf{A}$, at first glance it seems that because $\mathbf{A}^T=\mathbf{U}^T \mathbf{L}^T$, it might be that $\mathbf{U}$ and $\mathbf{L}$ are transposes of one another. But that's not possible in general, because we required only $\mathbf{L}$ to have ones on the diagonal, and that broke the symmetry. However, it's straightforward to restore the symmetry, as is demonstrated in {doc}`demos/structure-symm`.

If no pivoting is done, then the symmetric Gaussian elimination process yields $\mathbf{A}=\mathbf{L}\mathbf{D}\mathbf{L}^T$, where $\mathbf{L}$ is unit lower triangular and $\mathbf{D}$ is diagonal. In practice we don't actually have to carry out any arithmetic in the upper triangle as we work, since the operations are the mirror image of those in the lower triangle. As a result, it can be shown that LDL$^T$ factorization takes about half as much work as the standard LU, or $\sim \frac{1}{3}n^3$ flops.

In the general case we know that row pivoting is necessary to stabilize LU factorization. Pivoting is also needed to keep LDL$^T$ stable, but it has to be done symmetrically. We won't go into the details, as the resulting factorization isn't used very often. Instead, we'll explore the case when the matrix also possesses another important property.

(sec-SPD)=

## Symmetric positive definite matrices

Suppose that $\mathbf{A}$ is $n\times n$ and $\mathbf{x}\in\mathbb{R}^n$. Observe that $\mathbf{x}^T\mathbf{A}\mathbf{x}$ is the product of $1\times n$, $n\times n$, and $n\times 1$ matrices, so it is a scalar, sometimes referred to as a **quadratic form**.
In componentwise terms it becomes

```{math}
  :label: quadratic-form
  \mathbf{x}^T\mathbf{A}\mathbf{x} = \sum_{i=1}^n \sum_{j=1}^n A_{ij}x_ix_j.
```

```{index} symmetric positive definite matrix
```

```{index} see: SPD matrix; symmetric positive definite matrix
```

A real matrix $\mathbf{A}$ is called a {term}`symmetric positive definite matrix` (or SPD matrix) if it is symmetric and

```{math}
  :label: SPD-def
  \mathbf{x}^T \mathbf{A} \mathbf{x} > 0 \; \text{for all nonzero $\mathbf{x}\in\mathbb{R}^n$.}
```

The definiteness property is usually difficult to check directly from the definition. There are equivalent conditions, though; for instance, a symmetric matrix is positive definite if and only if its eigenvalues are all real positive numbers. SPD matrices have important properties and appear in applications in which the definiteness is known for theoretical reasons.

Let us consider what definiteness means to the LDL$^T$ factorization (itself the adaptation of LU factorization to symmetric matrices). We compute

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

```{index} matrix factorization; Cholesky
```

```{index} pivoting
```

Recall that the unpivoted LDL$^T$ (like unpivoted LU) factorization is not stable and not even always possible. However, in the SPD case one can prove that pivoting is not necessary for the existence nor the stability of the factorization $\mathbf{A}=\mathbf{R}^T\mathbf{R}$, which is known as the {term}`Cholesky factorization`. The elimination process is readily adapted into an algorithm for Cholesky factorization. Like LDL$^T$, the Cholesky factorization requires ${\sim\frac{1}{3}} n^3$ flops asymptotically in the $n\times n$ case, half as many as standard LU factorization.

```{proof:example} Julia demo
:class: demo
{doc}`demos/structure-cholesky`
```

The speed and stability of the Cholesky factorization make it the top choice for solving linear systems with SPD matrices. As a side benefit, the Cholesky algorithm fails, in the form of trying to take a negative square root or divide by zero, if and only if the matrix $\mathbf{A}$ is indefinite (i.e., not positive definite). In practice this is the best way to test the definiteness of a symmetric matrix about which nothing else is known.

## Exercises

1. ✍  For each matrix, use {eq}`diag-dominant` to determine whether it is diagonally dominant.

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
    \end{bmatrix}
    \quad \mathbf{C} =
    \begin{bmatrix}
    2  & -1 & 0  & 0      \\
    -1 & 2  & -1 & 0      \\
    0  & -1 & 2  & -1     \\
    0  & 0  & -1 & 2
    \end{bmatrix}.
    ```

2. ⌨ For each matrix, use inspection or `cholesky` in Julia to determine whether it is SPD.

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

    %    $A$ is the only SPD matrix. $B$ is not symmetric and $C$ causes
    %    "chol" to fail.
    %%
    %chol( [1 0 -1; 0 4 5; -1 5 10] )

    %%
    %chol([1 0 1; 0 4 5; -1 5 10])
    %%
    % This matrix is _not_ SPD, because it's not even symmetric, which |chol| doesn't check.

    %%
    %chol([1 0 1;0 4 5;1 5 1])

3. ✍ Show that the diagonal entries of a positive definite matrix are positive numbers. (Hint: Use special cases of {eq}`SPD-def`.)

4. ⌨ Using {ref}`function-lufact` as a guide, write a function

    ``` julia
    function luband(A,upper,lower)
    ```

    that accepts upper and lower bandwidth values and returns LU factors (without pivoting) in a way that avoids doing arithmetic using the locations that are known to stay zero.

5. ⌨ In this problem you will explore backslash and how it handles banded matrices. To do this you will generate tridiagonal matrices using the following code (supposing that `n` and `a` are defined):

    ``` julia
    A = triu( tril(rand(n,n),1), -1)
    A[1:n+1:end] .= a
    ```

    The result is $n\times n$, with each entry on the sub- and superdiagonals chosen randomly from $(0,1)$ and each diagonal entry equalling `a`.

    **(a)** Write a script that solves 200 linear systems whose matrices are generated as above, with $n=1000$ and $\mbox{"a"}=2$. (Use randomly generated right-hand side vectors.) Record the total time used by the solution process `A\b` only, using the built-in `@elapsed` timer. 

    **(b)** Repeat the experiment of part (a), but add the command `A=sparse(A)` right after the two lines above. How does the timing change?

    **(c)** Repeat parts (a) and (b) with $a=1$. In which case is there a major change?

    **(d)** Based on these observations, state a hypothesis on how backslash solves tridiagonal linear systems given in standard dense form and in sparse form. (Hint: What is the mathematical
    significance of $a=2$ versus $a=1$?)

6. ✍ A matrix $\mathbf{A}$ is called *skew-symmetric* if $\mathbf{A}^T=-\mathbf{A}$. Explain why unpivoted LU factorization of a skew-symmetric matrix is impossible.

    (problem-ATAisspd)=
7. ✍ Prove that if $\mathbf{A}$ is any real nonsingular square matrix, then $\mathbf{A}^T\mathbf{A}$ is SPD.

    %%
    % Let $x$ be any vector of the right length. Then
    %
    % $$x^T(A^TA)x = (Ax)^T(Ax) = \|Ax\|_2^2 \ge 0.$$
    %%
    % The only way to get equality is if $Ax=0$. But since $A$ is nonsingular, this means $x=0$.
