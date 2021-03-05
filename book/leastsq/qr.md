# The QR factorization

```{index} orthogonal; vectors
```

An important property of some groups of vectors is called **orthogonality**. We say that two vectors $\mathbf{u}$ and $\mathbf{v}$ in $\mathbb{R}^n$ are {term}`orthogonal` if $\mathbf{u}^T\mathbf{v}=0$. For $n=2$ or $n=3$ this means the vectors are perpendicular. We say that a collection of vectors $\mathbf{q}_1,\ldots,\mathbf{q}_k$ is orthogonal if

```{math}
:label: orthogonality
\mathbf{q}_i^T\mathbf{q}_j = 0 \quad \text{whenever $i\neq j$.}
```

```{index} orthonormal
```

If also $\mathbf{q}_i^T\mathbf{q}_i=1$ for all $i=1,\ldots,n$, we say the vectors are {term}`orthonormal`.

```{index} inner product
```

Orthogonal vectors are convenient theoretically and computationally. In theoretical calculations they make many terms of inner products vanish. For example, if $\mathbf{q}_1$ and $\mathbf{q}_2$ are orthogonal, then (in the 2-norm)

```{math}
:label: orthosubtract
\| \mathbf{q}_1 - \mathbf{q}_2 \|^2 = (\mathbf{q}_1-\mathbf{q}_2)^T(\mathbf{q}_1-\mathbf{q}_2)
= \mathbf{q}_1^T\mathbf{q}_1 - 2 \mathbf{q}_1^T\mathbf{q}_2 + \mathbf{q}_2^T\mathbf{q}_2
= \|\mathbf{q}_1\|^2 + \|\mathbf{q}_2\|^2.
```

This calculation is also the key to the computational attractiveness of orthogonality. {numref}`fig-nonorthogonal` shows how nonorthogonal vectors can allow a multidimensional version of subtractive cancellation, in which $\|\mathbf{x}-\mathbf{y}\|$ is much smaller than $\|\mathbf{x}\|$ and $\|\mathbf{y}\|$.

```{figure} figures/nonorthogonal.svg
:name: fig-nonorthogonal
Cancellation and orthogonality.
```

```{margin}
Addition and subtraction of vectors are guaranteed to be well conditioned when the vectors are orthogonal.
```

```{index} subtractive cancellation
```

As the figure and {eq}`orthosubtract` show, orthogonal vectors do not allow this phenomenon. In other words, the addition and subtraction of vectors are guaranteed to be well conditioned when the vectors are orthogonal.

## Orthogonal and ONC matrices

Statements about orthogonal vectors are often made more easily in matrix form. Let $\mathbf{Q}$ be an $n\times k$ matrix whose columns $\mathbf{q}_1$, \ldots, $\mathbf{q}_k$ are orthogonal vectors. The orthogonality conditions {eq}`orthogonality` become simply that $\mathbf{Q}^T\mathbf{Q}$ is a diagonal matrix, since

```{math}
\mathbf{Q}^T \mathbf{Q} =
\begin{bmatrix}
\mathbf{q}_1^T \\[1mm] \mathbf{q}_2^T \\ \vdots \\ \mathbf{q}_k^T
\end{bmatrix}
\begin{bmatrix}
\mathbf{q}_1 & \mathbf{q}_2 &  \cdots & \mathbf{q}_k
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{q}_1^T\mathbf{q}_1 & \mathbf{q}_1^T\mathbf{q}_2 & \cdots & \mathbf{q}_1^T\mathbf{q}_k \\[1mm]
\mathbf{q}_2^T\mathbf{q}_1 & \mathbf{q}_2^T\mathbf{q}_2 & \cdots & \mathbf{q}_2^T\mathbf{q}_k \\
\vdots & \vdots & & \vdots \\
\mathbf{q}_k^T\mathbf{q}_1 & \mathbf{q}_k^T\mathbf{q}_2 & \cdots & \mathbf{q}_k^T\mathbf{q}_k
\end{bmatrix}.
```

```{index} ONC matrix
```

If the columns of $\mathbf{Q}$ are orthonormal, then $\mathbf{Q}^T\mathbf{Q}$ is the $k\times k$ identity matrix. This is such an important property that we will break with common practice here and give this type of matrix a name: an {term}`ONC matrix` is one whose columns are an orthonormal set of vectors. We summarize their important properties here.

(theorem-ONC)=

````{prf:theorem} (ONC matrix)
Suppose $\mathbf{Q}$ is a real $n\times k$ ONC matrix (matrix with orthonormal columns). Then:

1. $\mathbf{Q}^T\mathbf{Q} = \mathbf{I}$ ($k\times k$ identity).
2. $\| \mathbf{Q}\mathbf{x} \|_2 = \| \mathbf{x} \|_2$ for all $k$-vectors $\mathbf{x}$.
3. $\| \mathbf{Q} \|_2=1$.
````

````{prf:proof}
The first part is derived above. The second part follows a pattern that has become well established by now:

```{math}
\| \mathbf{Q}\mathbf{x} \|_2^2 = (\mathbf{Q}\mathbf{x})^T(\mathbf{Q}\mathbf{x}) = \mathbf{x}^T \mathbf{Q}^T \mathbf{Q} \mathbf{x} = \mathbf{x}^T \mathbf{I} \mathbf{x} = \| \mathbf{x} \|_2^2.
```

The last part of the theorem is left to the exercises.
````

```{index} orthogonal; matrix
```

Of particular interest is a *square* ONC matrix, for which $\mathbf{Q}^T\mathbf{Q}=\mathbf{I}$, where all three matrices are $n\times n$. Hence $\mathbf{Q}^{-1}=\mathbf{Q}^T$. Such a matrix is called an {term}`orthogonal matrix`.[^ortho] These matrices have properties beyond the general ONC type. The proofs of these are left to the exercises.

[^ortho]: Confusingly, a square matrix whose columns are orthogonal is not necessarily an orthogonal matrix; the columns must be orthonormal, which is a stricter condition.

(theorem-orthogmatrix)=

````{prf:theorem} (Orthogonal matrix)
Suppose $\mathbf{Q}$ is an $n\times n$ real orthogonal matrix. Then:
1. $\mathbf{Q}^T$ is also an orthogonal matrix.
1. $\kappa(\mathbf{Q})=1$ in the 2-norm.
1. For any other $n\times n$ matrix $\mathbf{A}$, $\| \mathbf{A}\mathbf{Q} \|_2=\| \mathbf{A} \|_2$.
1. If $\mathbf{U}$ is another $n\times n$ orthogonal matrix, then $\mathbf{Q}\mathbf{U}$ is also orthogonal.
````

## Orthogonal factorization

```{margin}
The QR factorization plays a role in linear least squares analogous to the role of LU factorization in linear systems.
```

```{index} matrix factorization; QR
```

Now we come to another important way to factor a matrix, the {term}`QR factorization`. As we will show below, the QR factorization plays a role in linear least squares analogous to the role of LU factorization in linear systems.

(theorem-QR)=

````{prf:theorem}
  Every real $m\times n$ matrix $\mathbf{A}$ ($m\ge n$) can be written as $\mathbf{A}=\mathbf{Q}\mathbf{R}$, where $\mathbf{Q}$ is an $m\times m$ orthogonal matrix and $\mathbf{R}$ is an $m\times n$ upper triangular matrix.
````

In most introductory books on linear algebra, the QR factorization is derived through a process known as **Gram--Schmidt orthogonalization**. However, while it is an important tool for theoretical work, the Gram--Schmidt process is numerically unstable. We will consider an alternative construction in the next section.

When $m$ is much larger than $n$, which is often the case, there is a compressed form of the factorization that is more efficient. In the product

```{math}
\mathbf{A} =
\begin{bmatrix}
\mathbf{q}_1 & \mathbf{q}_2 & \cdots & \mathbf{q}_m
\end{bmatrix}
\begin{bmatrix}
r_{11} & r_{12} & \cdots & r_{1n} \\
0 & r_{22} & \cdots &  r_{2n} \\
\vdots & & \ddots & \vdots\\
0 & 0 & \cdots & r_{nn} \\
0 & 0 & \cdots & 0 \\
\vdots & \vdots &  & \vdots \\
0 & 0 & \cdots & 0
\end{bmatrix},
```

the vectors $\mathbf{q}_{n+1},\ldots,\mathbf{q}_m$ always get multiplied by zero. Nothing about $\mathbf{A}$ is lost if we delete them and reduce the factorization to the equivalent form

```{math}
:label: economyqr
\mathbf{A} =
\begin{bmatrix}
\mathbf{q}_1 & \mathbf{q}_2 & \cdots & \mathbf{q}_n
\end{bmatrix}
\begin{bmatrix}
r_{11} & r_{12} & \cdots & r_{1n} \\
0 & r_{22} & \cdots &  r_{2n} \\
\vdots & & \ddots & \vdots\\
0 & 0 & \cdots & r_{nn}
\end{bmatrix} = \hat{\mathbf{Q}} \hat{\mathbf{R}},
```

````{prf:example} Julia demo
:class: demo
{doc}`demos/qr-qrfact`
````

in which $\hat{\mathbf{Q}}$ is an $m\times n$ ONC matrix and $\hat{\mathbf{R}}$ is $n\times n$ and upper triangular. We refer to this as a **thin QR factorization**, as the number of columns in $\hat{\mathbf{Q}}$ is $n$ rather than $m$. Julia returns either type of QR factorization from the {term}`qr` command.

## Least squares and QR

```{index} linear least squares problem
```

```{index} normal equations
```

If we substitute the thin factorization {eq}`economyqr` into the normal equations {eq}`normaleqns`, we can simplify expressions a great deal.

```{math}
\begin{split}
  \mathbf{A}^T\mathbf{A} \mathbf{x} &= \mathbf{A}^T \mathbf{b} \\
  \hat{\mathbf{R}}^T \hat{\mathbf{Q}}^T \hat{\mathbf{Q}} \hat{\mathbf{R}} \mathbf{x} &= \hat{\mathbf{R}}^T \hat{\mathbf{Q}}^T \mathbf{b} \\
  \hat{\mathbf{R}}^T \hat{\mathbf{R}} \mathbf{x}& = \hat{\mathbf{R}}^T \hat{\mathbf{Q}}^T \mathbf{b}.
\end{split}
```

In order to have the normal equations be well posed, we require that $\mathbf{A}$ is not rank-deficient (as proved in a [theorem](theorem-ATA)). This is enough to guarantee that $\hat{\mathbf{R}}$ is nonsingular (see [this exercise](problem-nonsingR). Therefore, its transpose is nonsingular as well, and we arrive at

```{math}
:label: lsqr
\hat{\mathbf{R}} \mathbf{x}=\hat{\mathbf{Q}}^T \mathbf{b}.
```

```{margin}
The solution of least squares problems via QR factorization does not suffer from the instability seen when the normal equations are solved by Cholesky factorization.
```

This is a triangular $n\times n$ linear system that is easily solved by backward substitution, as demonstrated in {ref}`function-lsqrfact`. The function itself is superfluous, however, as this is essentially the algorithm used internally by Julia when \verb+A\b+ is called. Most importantly, even though we derived {eq}`lsqr` from the normal equations, the solution of least squares problems via QR factorization does not suffer from the instability seen when the normal equations are solved directly using Cholesky factorization.

(function-lsqrfact)=

```{proof:function} lsqrfact

```{code-block} julia
:lineno-start: 1
"""
lsqrfact(A,b)

Solve a linear least squares problem by QR factorization. Returns
the minimizer of ||b-Ax||.
"""
function lsqrfact(A,b)

Q,R = qr(A)
c = Q'*b
x = backsub(R,c)

return x
end
```

## Exercises

1. ✍ Prove part 3 of [the ONC matrix theorem](theorem-ONC).

2. ✍ Prove [the orthogonal matrix theorem](theorem-orthogmatrix). For the third part, use the definition of the 2-norm as an induced matrix norm, then apply some of our other results as needed.

3. ⌨ Let $t_0,\ldots,t_m$ be $m+1$ equally spaced points in $[-1,1]$. Let $\mathbf{A}$ be the matrix in {eq}`vandersystemrect` for $m=400$ and fitting by polynomials of degree less than 5. Find the thin QR factorization of $\mathbf{A}$, and, on a single graph, plot every column of $\hat{\mathbf{Q}}$ as a function of the vector $t$.

    <!-- t = linspace(-1,1,401)';
    A = [ t.^0 t.^1 t.^2 t.^3 t.^4 ];
    [Q,R] = qr(A,0);
    plot(t,Q) -->

    (problem-nonsingR)=
4. ✍ Prove that if the $m\times n$ ($m\ge n$) matrix $\mathbf{A}$ is not rank-deficient, then the factor $\hat{\mathbf{R}}$ of the thin QR factorization is nonsingular. (Hint: Suppose on the contrary that $\hat{\mathbf{R}}$ is singular. Show using the factored form of $\mathbf{A}$ that this would imply that $\mathbf{A}$ is rank-deficient.)

5. ✍ Let $\mathbf{A}$ be $m\times n$ with $m>n$. Show that if $\mathbf{A}=\mathbf{Q}\mathbf{R}$ is a QR factorization and $\mathbf{R}$ has rank $n$, then $\mathbf{A}^+=\mathbf{R}^+\mathbf{Q}^T$.

6. ✍ Let $\mathbf{A}$ be $m\times n$ with $m>n$. Show that if $\mathbf{A}=\hat{\mathbf{Q}}\hat{\mathbf{R}}$ is a QR factorization and $\hat{\mathbf{R}}$ is nonsingular, then $\mathbf{A}^+=\hat{\mathbf{R}}^{-1}\hat{\mathbf{Q}}^T$.

    %%
    % The QR factorization can be written as follows.
    %
    % $$A = QR = \left[ \hat{Q} \ Q_0 \right] \left[ \begin{array}{c} \hat{R} \\ 0 \end{array} \right]$$
    %

    %%
    % Note that $A$ is $m$-by-$n$, $\hat{Q}$ is $m$-by-$n$ and $\hat{R}_0$ is $n$-by-$n$; the others fill out
    % the dimensions with $Q_0$ being $m$-by-$(m-n)$ and the zero matrix is $(m-n)$-by-$n$.  Because of that
    % zero matrix, we can reconstruct $A$ with the compressed QR facorization: $A = \hat{Q}\hat{R}$.  Now use this
    % in the formula for the pseudoinverse.
    %
    % $$A^+  =  (A^T A)^{-1} A^T  =  [ (\hat{Q}\hat{R})^T \hat{Q}\hat{R} ]^{-1} [ \hat{Q}\hat{R} ]^T$$

    %%
    % $$=  [ \hat{R}^T \hat{Q}^T \hat{Q} \hat{R} ]^{-1} \hat{R}^T \hat{Q}^T$$

    %%
    % $$=  [ \hat{R}^T \hat{R} ]^{-1} \hat{R}^T \hat{Q}^T$$

    %%
    % $$= \hat{R}^{-1} [ \hat{R}^T ]^{-1} \hat{R}^T \hat{Q}^T$$
    %%
    % $$= \hat{R}^{-1} I \hat{Q}^T = \hat{R}^{-1} \hat{Q}^T.$$

7. ⌨ Repeat the [census fitting problem](problem-fitcensus), but use thin QR factorization rather than the backslash to solve the least-squares problem.

    (problem-orthoproj)=
8. ✍ The matrix $\mathbf{P}=\hat{\mathbf{Q}} \hat{\mathbf{Q}}^T$ derived from the thin QR factorization has some interesting and important properties.

    **(a)** Show that $\mathbf{P}=\mathbf{A}\mathbf{A}^+$.

    **(b)** Prove that $\mathbf{P}^2=\mathbf{P}$. (This property defines a *projection matrix*.)

    **(c)** Any vector $\mathbf{x}$ may be written as $\mathbf{x}=\mathbf{u}+\mathbf{v}$, where $\mathbf{u}=\mathbf{P}\mathbf{x}$ and $\mathbf{v}=(\mathbf{I}-\mathbf{P})\mathbf{x}$. Prove that for $\mathbf{P} =\hat{\mathbf{Q}} \hat{\mathbf{Q}}^T$, $\mathbf{u}$ and $\mathbf{v}$ are orthogonal. (Together with part (b), this means that $\mathbf{P}$ is an *orthogonal projector*.)
