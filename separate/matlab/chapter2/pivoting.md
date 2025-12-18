---
numbering:
  enumerator: 2.6.%s
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
---
```{code-cell}
:tags: [remove-cell]
clear all
format short
set(0, 'defaultaxesfontsize', 12)
set(0, 'defaultlinelinewidth', 1.5)
set(0, 'defaultFunctionLinelinewidth', 1.5)
set(0, 'defaultscattermarkerfacecolor', 'flat')
gcf;
set(gcf, 'Position', [0 0 600 350])
addpath FNC-matlab
```

(section-linsys-pivoting)=

# Row pivoting

```{index} matrix factorization; LU
```

As mentioned in {numref}`section-linsys-lu`, the $\mathbf{A}=\mathbf{L}\mathbf{U}$ factorization is not stable for every nonsingular $\mathbf{A}$. Indeed, the factorization does not always even exist.

::::{prf:example} Failure of naive LU factorization
:label: demo-pivoting-fail

Here is a previously encountered matrix that factors well.

```{code-cell}
A = [
    2 0 4 3
    -4 5 -7 -10
    1 15 2 -4.5
    -2 0 2 -13
    ];
[L, U] = lufact(A);
L
```

If we swap the second and fourth rows of $\mathbf{A}$, the result is still nonsingular. However, the factorization now fails.

```{code-cell}
A([2, 4], :) = A([4, 2], :);    % swap rows 2 and 4
[L, U] = lufact(A);
L
```

```{index} MATLAB; NaN
```

The presence of `NaN` in the result indicates that some impossible operation was required. The source of the problem is easy to locate. We can find the first outer product in the factorization just fine:

```{code-cell}
U(1, :) = A(1, :);
L(:, 1) = A(:, 1) / U(1, 1)
A = A - L(:, 1) * U(1, :)
```

The next step is `U(2, :) = A(2, :)`, which is also OK. But then we are supposed to divide by `U(2, 2)`, which is zero. The algorithm cannot continue.

::::

In {numref}`section-linsys-lu` we remarked that LU factorization is equivalent to Gaussian elimination with no row swaps. However, those swaps are necessary in situations like those encountered in @demo-pivoting-fail, in order to avoid division by zero. We will find a modification of the outer product procedure that allows us to do the same thing.

## Choosing a pivot

```{index} ! pivoting
```

The diagonal element of $\mathbf{U}$ that appears in the denominator of line 17 of {numref}`Function {number} <function-lufact>` is called the **pivot element** of its column. In order to avoid a zero pivot, we will use the largest available element in the column we are working on as the pivot. This technique is known as **row pivoting**.

```{prf:definition} Row pivoting
:label: definition-pivoting
When performing elimination in column $j$, choose as the pivot the element in column $j$ that is largest in absolute value. (In case of ties, choose the lowest row index.)
```

::::{prf:example} Row pivoting in LU factorization
:label: demo-pivoting-fix

Here is the trouble-making matrix from @demo-pivoting-fail.

```{code-cell}
A_1 = [2 0 4 3; -2 0 2 -13; 1 15 2 -4.5; -4 5 -7 -10]
```

```{index} ! MATLAB; max, ! MATLAB; \~
```

We now find the largest candidate pivot in the first column. We don't care about sign, so we take absolute values before finding the max.
```{tip}
:class: dropdown
The second output of `max` returns the location of the largest element of a vector. The `~` symbol is used to ignore the value of the first output.
```


```{code-cell}
[~, i] = max( abs(A_1(:, 1)) ) 
```

This is the row of the matrix that we extract to put into $\mathbf{U}$. That guarantees that the division used to find $\boldsymbol{\ell}_1$ will be valid.

```{code-cell}
L = zeros(4, 4);
U = zeros(4, 4);
U(1, :) = A_1(i, :);
L(:, 1) = A_1(:, 1) / U(1, 1);
A_2 = A_1 - L(:, 1) * U(1, :)
```

Observe that $\mathbf{A}_2$ has a new zero row and zero column, but the zero row is the fourth rather than the first. However, we forge on by using the largest possible pivot in column 2 for the next outer product.

```{code-cell}
[~, i] = max( abs(A_2(:, 2)) )
U(2, :) = A_2(i, :);
L(:, 2) = A_2(:, 2) / U(2, 2);
A_3 = A_2 - L(:, 2) * U(2, :)
```

Now we have zeroed out the third row as well as the second column. We can finish out the procedure.

```{code-cell}
[~, i] = max( abs(A_3(:, 3)) ) 
U(3, :) = A_3(i, :);
L(:, 3) = A_3(:, 3) / U(3, 3);
A_4 = A_3 - L(:, 3) * U(3, :)
```

```{code-cell}
[~, i] = max( abs(A_4(:, 4)) ) 
U(4, :) = A_4(i, :);
L(:, 4) = A_4(:, 4) / U(4, 4);
```

We do have a factorization of the original matrix:

```{code-cell}
A_1 - L * U
```

And $\mathbf{U}$ has the required structure:

```{code-cell}
U
```

However, the triangularity of $\mathbf{L}$ has been broken.

```{code-cell}
L
```

::::

We will return to the loss of triangularity in $\mathbf{L}$ momentarily. First, though, there is a question left to answer: what if at some stage, all the elements of the targeted column are zero, i.e., there are no available pivots? Fortunately that loose end ties up nicely, although a proof is a bit beyond our scope here.

```{prf:theorem} Row pivoting
:label: theorem-pivot
The row-pivoted LU factorization runs to completion if and only if the original matrix is invertible.
```

A linear system with a singular matrix has either no solution or infinitely many solutions. Either way, a technique other than LU factorization is needed to handle it.

## Permutations

Even though the resulting $\mathbf{L}$ in @demo-pivoting-fix is no longer of unit lower triangular form, it is close. In fact, all that is needed is to reverse the order of its rows.

::::{prf:example} Pivoting as row permutation
:label: demo-pivoting-permute

Here again is the matrix from @demo-pivoting-fix.

```{code-cell}
A = [2 0 4 3; -2 0 2 -13; 1 15 2 -4.5; -4 5 -7 -10]
```

As the factorization proceeded, the pivots were selected from rows 4, 3, 2, and finally 1. If we were to put the rows of $\mathbf{A}$ into that order, then the algorithm would run exactly like the plain LU factorization from {numref}`section-linsys-lu`. 

```{code-cell}
B = A([4, 3, 2, 1], :);
[L, U] = lufact(B);
```

We obtain the same $\mathbf{U}$ as before:

```{code-cell}
U
```

And $\mathbf{L}$ has the same rows as before, but arranged into triangular order:

```{code-cell}
L
```

::::

In principle, if the permutation of rows implied by the pivot locations is applied all at once to the original $\mathbf{A}$, no further pivoting is needed. In practice, this permutation cannot be determined immediately from the original $\mathbf{A}$; the only way to find it is to run the algorithm. Having obtained it at the end, though, we can use it to state a simple relationship.

```{index} ! matrix factorization; pivoted LU, ! PLU factorization
```

::::{prf:definition} PLU factorization
:label: definition-plu
Given $n\times n$ matrix $\mathbf{A}$, the {term}`PLU factorization` is a unit lower triangular $\mathbf{L}$, an upper triangular $\mathbf{U}$, and a permutation $i_1,\ldots,i_n$ of the integers $1,\ldots,n$, such that

$$\tilde{\mathbf{A}} = \mathbf{L}\mathbf{U},$$

where rows $1,\ldots,n$ of $\tilde{\mathbf{A}}$ are rows $i_1,\ldots,i_n$ of $\mathbf{A}$.
::::

{numref}`Function {number} <function-plufact>` shows our implementation of PLU factorization.[^PLU]

[^PLU]: Because unpivoted LU factorization is not useful, in practice the term *LU factorization* mostly refers to pivoted LU.

``````{prf:algorithm} plufact
:label: function-plufact

```{literalinclude} ../FNC_matlab/plufact.m
:linenos: true
:language: matlab
```
```
``````

Ideally, the PLU factorization takes $\sim \frac{2}{3}n^3$ flops asymptotically, just like LU without pivoting. The implementation in {numref}`Function {number} <function-plufact>` does not achieve this optimal flop count, however. Like {numref}`Function {number} <function-lufact>`, it does unnecessary operations on structurally known zeros for the sake of being easier to understand.

## Linear systems

The output of {numref}`Function {number} <function-plufact>` is a factorization of a row-permuted $\mathbf{A}$. Therefore, given a linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$, we have to permute $\mathbf{b}$ the same way before applying forward and backward substitution. This is equivalent to changing the order of the equations in a linear system, which does not affect its solution.

::::{prf:example} PLU factorization for solving linear systems
:label: demo-pivoting-usage

The third output of `plufact` is the permutation vector we need to apply to $\mathbf{A}$.

```{code-cell}
A = randi(20, 4, 4);
[L, U, p] = plufact(A);
A(p, :) - L * U    % should be ≈ 0
```

Given a vector $\mathbf{b}$, we solve $\mathbf{A}\mathbf{x}=\mathbf{b}$ by first permuting the entries of $\mathbf{b}$ and then proceeding as before.

```{code-cell}
b = rand(4, 1);
z = forwardsub(L, b(p));
x = backsub(U, z)
```

A residual check is successful:

```{code-cell}
b - A*x
```

::::

While {numref}`Function {number} <function-plufact>` is a serviceable implementation, it is not the gold standard for solving linear systems in practice. Each language offers its own built-in function for PLU factorization.

::::{prf:example} Built-in PLU factorization
:label: demo-pivoting-builtin

With the syntax `A \ b`, the matrix `A` is PLU-factored, followed by two triangular solves.

```{code-cell}
A = randn(500, 500);    % 500x500 with normal random entries
tic; for k=1:50; A \ rand(500, 1); end; toc
```

In {numref}`section-linsys-efficiency` we showed that the factorization is by far the most costly part of the solution process. A factorization object allows us to do that costly step only once per unique matrix. 

```{code-cell}
[L, U, p] = lu(A, 'vector');    % keep factorization result
tic
for k=1:50
    b = rand(500, 1);
    U \ (L \ b(p));
end
toc
```

::::

## Stability

```{index} stability
```

There is one detail of the row pivoting algorithm that might seem arbitrary: why choose the pivot of largest magnitude within a column, rather than, say, the uppermost nonzero in the column? The answer is numerical stability.

::::{prf:example} Stability of PLU factorization
:label: demo-pivoting-stable
Let

```{math}
:label: plu-stab-A
\mathbf{A} =
  \begin{bmatrix}
    -\epsilon & 1 \\ 1 & -1
  \end{bmatrix}.
```

If $\epsilon=0$, LU factorization without pivoting fails for $\mathbf{A}$. But if $\epsilon\neq 0$, we can go without pivoting, at least in principle.

We construct a linear system for this matrix with $\epsilon=10^{-12}$ and exact solution $[1, 1]$:

```{code-cell}
ep = 1e-12
A = [-ep 1; 1 -1];
b = A * [1; 1];
```

We can factor the matrix without pivoting and solve for $\mathbf{x}$.

```{code-cell}
[L, U] = lufact(A);
x = backsub( U, forwardsub(L, b) )
```

Note that we have obtained only about 5 accurate digits for $x_1$. We could make the result even more inaccurate by making $\epsilon$ even smaller:

```{code-cell}
ep = 1e-20; A = [-ep 1; 1 -1];
b = A * [1; 1];
[L, U] = lufact(A);
x = backsub( U, forwardsub(L, b) )
```

This effect is not due to ill conditioning of the problem—a solution with PLU factorization works perfectly:

```{code-cell}
A \ b
```

::::

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

``````{exercise}
:label: problem-pivoting-byhand

✍ Perform by hand the pivoted LU factorization of each matrix.

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
``````

``````{exercise}
:label: problem-pivoting-grouping

✍ Let $\mathbf{A}$ be a square matrix and $\mathbf{b}$ be a column vector of compatible length. 

Here is correct code to solve $\mathbf{A}\mathbf{x}=\mathbf{b}$:

``` matlab
[L, U, p] = lu(A, 'vector');
x = U \ (L \ b(p));
```

Suppose instead you replace the last line above with

``` matlab
x = U \ L \ b(p);
```

Mathematically in terms of $\mathbf{L}$, $\mathbf{U}$, $\mathbf{p}$, and $\mathbf{b}$, what vector is found?
``````

``````{exercise}
:label: problem-pivoting-flip

✍ Suppose that `A` is a $4\times 6$ matrix and you define

``` matlab
B = A(end:-1:1, end:-1:1)
```

Show that $\mathbf{B} = \mathbf{P} \mathbf{A} \mathbf{Q}$ for certain matrices $\mathbf{P}$ and $\mathbf{Q}$.
```````

``````{exercise}
:label: problem-pivoting-perminverse
✍ An $n\times n$ *permutation matrix* $\mathbf{P}$ is a reordering of the rows of an identity matrix such that $\mathbf{P} \mathbf{A}$  has the effect of moving rows $1,2,\ldots,n$ of $\mathbf{A}$ to new positions $i_1,i_2,\ldots,i_n$. Then $\mathbf{P}$ can be expressed as

```{math}
:numbered: false
\mathbf{P} = \mathbf{e}_{i_1}\mathbf{e}_1^T + \mathbf{e}_{i_2}\mathbf{e}_2^T + \cdots + \mathbf{e}_{i_n}\mathbf{e}_n^T.
```

**(a)** For the case $n=4$ and $i_1=3$, $i_2=2$, $i_3=4$, $i_4=1$, write out separately, as matrices, all four of the terms in the sum. Then add them together to find $\mathbf{P}$.

**(b)** Use the formula in the general case to show that $\mathbf{P}^{-1}=\mathbf{P}^T$. (Hint: To show that $\mathbf{Z}$ is the inverse of $\mathbf{P}$, show that $\mathbf{P}\mathbf{Z}=\mathbf{I}$ or $\mathbf{Z}\mathbf{P}=\mathbf{I}$.)
``````
