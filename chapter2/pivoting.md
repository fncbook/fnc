---
jupytext:
  cell_metadata_filter: -all
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: Julia 1.7.1
  language: julia
  name: julia-fast
---
```{code-cell}
:tags: [remove-cell]
using FundamentalsNumericalComputation
FNC.init_format()
```

(section-linsys-pivoting)=
# Row pivoting

```{index} matrix factorization; LU 
```

As mentioned in {numref}`section-linsys-lu`, the $\mathbf{A}=\mathbf{L}\mathbf{U}$ factorization is not stable for every nonsingular $\mathbf{A}$. Indeed, the factorization does not always even exist.

(demo-pivoting-fail)=
```{proof:demo}
```

```{raw} html
<div class='demo'>
```

```{raw} latex
%%start demo%%
```

Here is a previously encountered matrix, which factors well.

```{code-cell}
A = [2 0 4 3 ; -4 5 -7 -10 ; 1 15 2 -4.5 ; -2 0 2 -13];
L,U = FNC.lufact(A)
L
```

If we swap the second and fourth rows of $\mathbf{A}$, the result is still nonsingular. However, the factorization now fails.

```{code-cell}
A[[2,4],:] = A[[4,2],:]  
L,U = FNC.lufact(A)
L
```

```{index} Julia; NaN
```

The presence of `NaN` in the result indicates that some impossible operation was required. The source of the problem is easy to locate. We can find the first outer product in the factorization just fine:

```{code-cell}
U[1,:] = A[1,:]
L[:,1] = A[:,1]/U[1,1]
A -= L[:,1]*U[1,:]'
```

The next step is `U[2,:]=A[2,:]`, which is also OK. But then we are supposed to divide by `U[2,2]`, which is zero. The algorithm cannot continue.

```{raw} html
</div>
```

```{raw} latex
%%end demo%%
```

In {numref}`section-linsys-lu` we remarked that LU factorization is equivalent to Gaussian elimination with no row swaps. However, those swaps are necessary in situations like those encountered in {numref}`Demo {number} <demo-pivoting-fail>`, in order to avoid division by zero. We will find a modification of the outer product procedure that allows us to do the same thing.

## Choosing a pivot

```{index} ! pivoting
```

The diagonal element of $\mathbf{U}$ that appears in the denominator of line 17 of {numref}`Function {number} <function-lufact>` is called the **pivot element** of its column. In order to avoid a zero pivot, we will use the largest available element in the column we are working on as the pivot. This technique is known as **row pivoting**.

(rule-pivoting)=
```{proof:algorithm} Row pivoting
When performing elimination in column $j$, choose as the pivot the element in column $j$ that is largest in absolute value. (In case of ties, choose the lowest row index.)
```

(demo-pivoting-fix)=
```{proof:demo}
```

```{raw} html
<div class='demo'>
```

```{raw} latex
%%start demo%%
```

Here is the trouble-making matrix from {numref}`Demo {number} <demo-pivoting-fail>`.

```{code-cell}
A₁ = [2 0 4 3 ; -2 0 2 -13; 1 15 2 -4.5 ; -4 5 -7 -10]
```

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item}
:columns: 7

```{raw} latex
\begin{minipage}[t]{0.5\textwidth}
```
We now find the largest candidate pivot in the first column. We don't care about sign, so we take absolute values before finding the max.

```{raw} latex
\end{minipage}\hfill
```
:::
:::{grid-item-card}
:columns: 5

```{raw} latex
\begin{minipage}[t]{0.4\textwidth}\begin{mdframed}[default]\small
```
The `argmax` function returns the location of the largest element of a vector or matrix.
```{raw} latex
\end{mdframed}\end{minipage}
```
:::
::::


```{code-cell}
i = argmax( abs.(A₁[:,1]) ) 
```

This is the row of the matrix that we extract to put into $\mathbf{U}$. That guarantees that the division used to find $\boldsymbol{\ell}_1$ will be valid.

```{code-cell}
L,U = zeros(4,4),zeros(4,4)
U[1,:] = A₁[i,:]
L[:,1] = A₁[:,1]/U[1,1]
A₂ = A₁ - L[:,1]*U[1,:]'
```

Observe that $\mathbf{A}_2$ has a new zero row and zero column, but the zero row is the fourth rather than the first. However, we forge on by using the largest possible pivot in column 2 for the next outer product.

```{code-cell}
@show i = argmax( abs.(A₂[:,2]) ) 
U[2,:] = A₂[i,:]
L[:,2] = A₂[:,2]/U[2,2]
A₃ = A₂ - L[:,2]*U[2,:]'
```

Now we have zeroed out the third row as well as the second column. We can finish out the procedure.

```{code-cell}
@show i = argmax( abs.(A₃[:,3]) ) 
U[3,:] = A₃[i,:]
L[:,3] = A₃[:,3]/U[3,3]
A₄ = A₃ - L[:,3]*U[3,:]'
```

```{code-cell}
@show i = argmax( abs.(A₄[:,4]) ) 
U[4,:] = A₄[i,:]
L[:,4] = A₄[:,4]/U[4,4];
```

We do have a factorization of the original matrix:

```{code-cell}
A₁ - L*U
```

And $\mathbf{U}$ has the required structure:

```{code-cell}
U
```

However, the triangularity of $\mathbf{L}$ has been broken.

```{code-cell}
L
```

```{raw} html
</div>
```

```{raw} latex
%%end demo%%
```

We will return to the loss of triangularity in $\mathbf{L}$ momentarily. First, though, there is a question left to answer: what if at some stage, all the elements of the targeted column are zero, i.e., there are no available pivots? Fortunately that loose end ties up nicely, although a proof is a bit beyond our scope here.

(theorem-pivot)=
```{proof:theorem} Row pivoting
The row-pivoted LU factorization runs to completion if and only if the original matrix is invertible.
```

A linear system with a singular matrix has either no solution or infinitely many solutions. Either way, a technique other than LU factorization is needed to handle it.

## Permutations

Even though the resulting $\mathbf{L}$ in {numref}`Demo {number} <demo-pivoting-fix>` is no longer of unit lower triangular form, it is close. In fact, all that is needed is to reverse the order of its rows. 

(demo-pivoting-permute)=
```{proof:demo}
```

```{raw} html
<div class='demo'>
```

```{raw} latex
%%start demo%%
```

Here again is the matrix from {numref}`Demo {number} <demo-pivoting-fix>`.

```{code-cell}
A = [2 0 4 3 ; -2 0 2 -13; 1 15 2 -4.5 ; -4 5 -7 -10]
```

As the factorization proceeded, the pivots were selected from rows 4, 3, 2, and finally 1. If we were to put the rows of $\mathbf{A}$ into that order, then the algorithm would run exactly like the plain LU factorization from {numref}`section-linsys-lu`. 

```{code-cell}
B = A[ [4,3,2,1], : ]
L,U = FNC.lufact(B);
```

We obtain the same $\mathbf{U}$ as before:

```{code-cell}
U
```

And $\mathbf{L}$ has the same rows as before, but arranged into triangular order:

```{code-cell}
L
```

```{raw} html
</div>
```

```{raw} latex
%%end demo%%
```

In principle, if the permutation of rows implied by the pivot locations is applied all at once to the original $\mathbf{A}$, no further pivoting is needed. In practice, this permutation cannot be determined immediately from the original $\mathbf{A}$; the only way to find it is to run the algorithm. Having obtained it at the end, though, we can use it to state a simple relationship. 

```{index} ! matrix factorization; pivoted LU, ! PLU factorization
```

::::{proof:definition} PLU factorization
Given $n\times n$ matrix $\mathbf{A}$, the **PLU factorization** is a unit lower triangular $\mathbf{L}$, an upper triangular $\mathbf{U}$, and a permutation $i_1,\ldots,i_n$ of the integers $1,\ldots,n$, such that

$$\tilde{\mathbf{A}} = \mathbf{L}\mathbf{U},$$

where rows $1,\ldots,n$ of $\tilde{\mathbf{A}}$ are rows $i_1,\ldots,i_n$ of $\mathbf{A}$.
::::

{numref}`Function {number} <function-plufact>` shows our implementation of PLU factorization.[^PLU]

[^PLU]: Because unpivoted LU factorization is not useful, in practice the term *LU factorization* mostly refers to pivoted LU.

(function-plufact)=
````{proof:function} plufact
**Row=pivoted LU factorization**
```{code-block} julia
:lineno-start: 1
"""
    plufact(A)

Compute the PLU factorization of square matrix `A`, returning the
triangular factors and a row permutation vector.
"""
function plufact(A)
    n = size(A,1)
    L = zeros(n,n)
    U = zeros(n,n)
    p = fill(0,n)
    Aₖ = float(copy(A))

    # Reduction by outer products
    for k in 1:n-1
        p[k] = argmax(abs.(Aₖ[:,k]))
        U[k,:] = Aₖ[p[k],:]
        L[:,k] = Aₖ[:,k]/U[k,k]
        Aₖ -= L[:,k]*U[k,:]'
    end
    p[n] = argmax(abs.(Aₖ[:,n]))
    U[n,n] = Aₖ[p[n],n]
    L[:,n] = Aₖ[:,n]/U[n,n]
    return LowerTriangular(L[p,:]),U,p
end
```
````

Ideally, the PLU factorization takes $\sim \frac{2}{3}n^3$ flops asymptotically, just like LU without pivoting. The implementation in {numref}`Function {number} <function-plufact>` does not achieve this optimal flop count, however. Like {numref}`Function {number} <function-lufact>`, it does unnecessary operations on structurally known zeros for the sake of being easier to understand.
## Linear systems

The output of {numref}`Function {number} <function-plufact>` is a factorization of a row-permuted $\mathbf{A}$. Therefore, given a linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$, we have to permute $\mathbf{b}$ the same way before applying forward and backward substitution. This is equivalent to changing the order of the equations in a linear system, which does not affect its solution.

(demo-pivoting-usage)=
:::{proof:demo}
:::

```{raw} html
<div class='demo'>
```

```{raw} latex
%%start demo%%
```

The third output of `plufact` is the permutation vector we need to apply to $\mathbf{A}$.

```{code-cell}
A = rand(1:20,4,4)
L,U,p = FNC.plufact(A)
A[p,:] - L*U   # should be ≈ 0
```

Given a vector $\mathbf{b}$, we solve $\mathbf{A}\mathbf{x}=\mathbf{b}$ by first permuting the entries of $\mathbf{b}$ and then proceeding as before.

```{code-cell}
b = rand(4)
z = FNC.forwardsub(L,b[p])
x = FNC.backsub(U,z)
```

A residual check is successful:

```{code-cell}
b - A*x
```

```{raw} html
</div>
```

```{raw} latex
%%end demo%%
```

The `lu` function from the built-in package `LinearAlgebra` returns the same three outputs as {numref}`Function {number} <function-plufact>`. If you only request one output, it will be a factorization object that can be used with a backslash. This is useful when you want to solve with multiple versions of $\mathbf{b}$ but do the factorization only once.

(demo-pivoting-builtin)=
:::{proof:demo}
:::

```{raw} html
<div class='demo'>
```

```{raw} latex
%%start demo%%
```

With the syntax `A\b`, the matrix `A` is PLU-factored, followed by two triangular solves.

```{code-cell}
A = randn(500,500)   # 500x500 with normal random entries
A\rand(500)          # force compilation
@elapsed for k=1:50; A\rand(500); end
```

In {numref}`section-linsys-efficiency` we showed that the factorization is by far the most costly part of the solution process. A factorization object allows us to do that costly step only once per unique matrix. 

```{code-cell}
factored = lu(A)     # store factorization result
factored\rand(500)   # force compilation
@elapsed for k=1:50; factored\rand(500); end
```

```{raw} html
</div>
```

```{raw} latex
%%end demo%%
```

## Stability

```{index} stability
```

There is one detail of the row pivoting algorithm that might seem arbitrary: why choose the pivot of largest magnitude in a column, rather than, say, the uppermost nonzero in the column? The answer is numerical stability.

(demo-pivoting-stable)=
:::{proof:demo}
:::

```{raw} html
<div class='demo'>
```

```{raw} latex
%%start demo%%
```

Let

```{math}
:label: plu-stab-A
\mathbf{A} =
  \begin{bmatrix}
    -\epsilon & 1 \\ 1 & -1
  \end{bmatrix}.
```

If $\epsilon=0$, LU factorization without pivoting fails for $\mathbf{A}$. But if $\epsilon\neq 0$, we can go without pivoting, at least in principle.

We construct a linear system for this matrix with $\epsilon=10^{-12}$ and exact solution $[1,1]$:

```{code-cell}
ϵ = 1e-12
A = [-ϵ 1;1 -1]
b = A*[1,1]
```

We can factor the matrix without pivoting and solve for $\mathbf{x}$.

```{code-cell}
L,U = FNC.lufact(A)
x = FNC.backsub( U, FNC.forwardsub(L,b) )
```

Note that we have obtained only about 5 accurate digits for $x_1$. We could make the result even more inaccurate by making $\epsilon$ even smaller:

```{code-cell}
ϵ = 1e-20; A = [-ϵ 1;1 -1]
b = A*[1,1]
L,U = FNC.lufact(A)
x = FNC.backsub( U, FNC.forwardsub(L,b) )
```

This effect is not due to ill conditioning of the problem—a solution with PLU factorization works perfectly:

```{code-cell}
A\b
```

```{raw} html
</div>
```

```{raw} latex
%%end demo%%
```

The factors of this $\mathbf{A}$ without pivoting are found to be

```{math}
:label: plu-stab-LU
  \mathbf{L} = 
  \begin{bmatrix}
    1 & 0 \\ -\epsilon^{-1} & 1 
  \end{bmatrix}, \qquad 
  \mathbf{U} = 
  \begin{bmatrix}
    -\epsilon & 1 \\ 0 & \epsilon^{-1}-1 
  \end{bmatrix}.
```

For reasons we will quantify in {numref}`section-linsys-condition-number`, the solution of $\mathbf{A}\mathbf{x}=\mathbf{b}$ is well-conditioned, but the problems of solving $\mathbf{L}\mathbf{z}=\mathbf{b}$ and $\mathbf{U}\mathbf{x}=\mathbf{z}$ have condition numbers essentially $1/\epsilon^2$ each. Thus, for small $\epsilon$, solution of the original linear system by unpivoted LU factorization is highly unstable. 

Somewhat surprisingly, solving $\mathbf{A}\mathbf{x}=\mathbf{b}$ via PLU factorization is technically also unstable. In fact, examples of unstable solutions are well-known, but they have been nonexistent in practice. While there is a lot of evidence and some reasoning about why this is the case, the situation is not completely understood. Yet PLU factorization remains the algorithm of choice for general linear systems.

## Exercises

1. ✍ Perform by hand the pivoted LU factorization of each matrix.
    
    **(a)** $\quad \displaystyle \begin{bmatrix}
    2 & 3 & 4 \\
    4 & 5 & 10 \\
    4 & 8 & 2
    \end{bmatrix},\qquad$
    **(b)** $\quad \displaystyle \begin{bmatrix}
    1 & 4 & 5 & -5 \\
    -1 & 0 & -1 & -5 \\
    1 & 3 & -1 & 2 \\
    1 & -1 & 5 & -1 
    \end{bmatrix}$.

2. ✍ Let $\mathbf{A}$ be a square matrix and $\mathbf{b}$ be a column vector of compatible length. Here is correct Julia code to solve $\mathbf{A}\mathbf{x}=\mathbf{b}$:

    ``` julia
    L,U,p = lu(A)
    x = U \ (L\b[p])
    ```

    Suppose instead you replace the last line above with

    ``` julia
    x = U \ L \ b[p]
    ```

    Mathematically in terms of $\mathbf{L}$, $\mathbf{U}$, $\mathbf{p}$, and $\mathbf{b}$, what vector is found?

3. ✍ Suppose that `A` is a $4\times 6$ matrix in Julia and you define

    ``` julia
    B = A[end:-1:1,end:-1:1]
    ```

    Show that $\mathbf{B} = \mathbf{P} \mathbf{A} \mathbf{Q}$ for certain matrices $\mathbf{P}$ and $\mathbf{Q}$.

    (problem-perminverse)=
4. ✍ An $n\times n$ *permutation matrix* $\mathbf{P}$ is a reordering of the rows of an identity matrix such that $\mathbf{P} \mathbf{A}$  has the effect of moving rows $1,2,\ldots,n$ of $\mathbf{A}$ to new positions $i_1,i_2,\ldots,i_n$. Then $\mathbf{P}$ can be expressed as
  
    ```{math}
    \mathbf{P} = \mathbf{e}_{i_1}\mathbf{e}_1^T + \mathbf{e}_{i_2}\mathbf{e}_2^T + \cdots + \mathbf{e}_{i_n}\mathbf{e}_n^T.
    ```

    **(a)** For the case $n=4$ and $i_1=3$, $i_2=2$, $i_3=4$, $i_4=1$, write out separately, as matrices, all four of the terms in the sum. Then add them together to find $\mathbf{P}$.

    **(b)** Use the formula in the general case to show that $\mathbf{P}^{-1}=\mathbf{P}^T$.
  

