---
numbering:
  enumerator: 2.6.%s
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---
```{code-cell}
:tags: [remove-cell]
from numpy import *
from scipy import linalg
from scipy.linalg import norm
from matplotlib.pyplot import *
from prettytable import PrettyTable
from timeit import default_timer as timer
import sys
sys.path.append('fncbook/')
import fncbook as FNC

# This (optional) block is for improving the display of plots.
# from IPython.display import set_matplotlib_formats
# set_matplotlib_formats("svg","pdf")
# %config InlineBackend.figure_format = 'svg'
rcParams["figure.figsize"] = [7, 4]
rcParams["lines.linewidth"] = 2
rcParams["lines.markersize"] = 4
rcParams['animation.html'] = "jshtml"  # or try "html5"
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
A = array([
    [2, 0, 4, 3],
    [-4, 5, -7, -10],
    [1, 15, 2, -4.5],
    [-2, 0, 2, -13]
    ])
L, U = FNC.lufact(A)
print(L)
```

If we swap the second and fourth rows of $\mathbf{A}$, the result is still nonsingular. However, the factorization now fails.

```{code-cell} 
A[[1, 3], :] = A[[3, 1], :]  
L, U = FNC.lufact(A)
print(L)
```

```{index} Python; NaN
```

The presence of `NaN` in the result indicates that some impossible operation was required. The source of the problem is easy to locate. We can find the first outer product in the factorization just fine:

```{code-cell}
U[0, :] = A[0, :]
L[:, 0] = A[:, 0] / U[0, 0]
A -= outer(L[:, 0],  U[0, :])
print(A)
```

The next step is `U[1, :] = A[1, :]`, which is also OK. But then we are supposed to divide by `U[1, 1]`, which is zero. The algorithm cannot continue.

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
A_1 = array([
    [2, 0, 4, 3],
    [-2, 0, 2, -13],
    [1, 15, 2, -4.5],
    [-4, 5, -7, -10]
    ])
```

We now find the largest candidate pivot in the first column. We don't care about sign, so we take absolute values before finding the max.
```{tip}
:class: dropdown
The `argmax` function returns the location of the largest element of a vector or matrix.
```


```{code-cell}
i = argmax( abs(A_1[:, 0]) )
print(i)
```

This is the row of the matrix that we extract to put into $\mathbf{U}$. That guarantees that the division used to find $\boldsymbol{\ell}_1$ will be valid.

```{code-cell}
L, U = eye(4), zeros((4, 4))
U[0, :] = A_1[i, :]
L[:, 0] = A_1[:, 0] / U[0, 0]
A_2 = A_1 - outer(L[:, 0], U[0, :])
print(A_2)
```

Observe that $\mathbf{A}_2$ has a new zero row and zero column, but the zero row is the fourth rather than the first. However, we forge on by using the largest possible pivot in column 2 for the next outer product.

```{code-cell}
i = argmax( abs(A_2[:, 1]) ) 
print(f"new pivot row is {i}")
U[1, :] = A_2[i, :]
L[:, 1] = A_2[:, 1] / U[1, 1]
A_3 = A_2 - outer(L[:, 1], U[1, :])
print(A_3)
```
Now we have zeroed out the third row as well as the second column. We can finish out the procedure.

```{code-cell}
i = argmax( abs(A_3[:, 2]) ) 
print(f"new pivot row is {i}")
U[2, :] = A_3[i, :]
L[:, 2] = A_3[:, 2] / U[2, 2]
A_4 = A_3 - outer(L[:, 2], U[2, :])
print(A_4)
```

```{code-cell}
i = argmax( abs(A_4[:, 3]) ) 
print(f"new pivot row is {i}")
U[3, :] = A_4[i, :]
L[:, 3] = A_4[:, 3] / U[3, 3];
```

We do have a factorization of the original matrix:

```{code-cell}
A_1 - L @ U
```

And $\mathbf{U}$ has the required structure:

```{code-cell}
print(U)
```

However, the triangularity of $\mathbf{L}$ has been broken.

```{code-cell}
print(L)
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
A = array([
    [2, 0, 4, 3],
    [-2, 0, 2, -13],
    [1, 15, 2, -4.5],
    [-4, 5, -7, -10]
    ])
```

As the factorization proceeded, the pivots were selected from rows 4, 3, 2, and finally 1 (with NumPy indices being one less). If we were to put the rows of $\mathbf{A}$ into that order, then the algorithm would run exactly like the plain LU factorization from {numref}`section-linsys-lu`. 

```{code-cell}
B = A[[3, 2, 1, 0], :]
L, U = FNC.lufact(B);
```

We obtain the same $\mathbf{U}$ as before:

```{code-cell}
print(U)
```

And $\mathbf{L}$ has the same rows as before, but arranged into triangular order:

```{code-cell}
print(L)
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

```{literalinclude} chapter02.py
:filename: plufact.py
:start-line: 51
:end-line: 71
:linenos: true
:language: python
```
``````

Ideally, the PLU factorization takes $\sim \frac{2}{3}n^3$ flops asymptotically, just like LU without pivoting. The implementation in {numref}`Function {number} <function-plufact>` does not achieve this optimal flop count, however. Like {numref}`Function {number} <function-lufact>`, it does unnecessary operations on structurally known zeros for the sake of being easier to understand.

## Linear systems

The output of {numref}`Function {number} <function-plufact>` is a factorization of a row-permuted $\mathbf{A}$. Therefore, given a linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$, we have to permute $\mathbf{b}$ the same way before applying forward and backward substitution. This is equivalent to changing the order of the equations in a linear system, which does not affect its solution.

::::{prf:example} PLU factorization for solving linear systems
:label: demo-pivoting-usage

The third output of `plufact` is the permutation vector we need to apply to $\mathbf{A}$.

```{code-cell}
A = random.randn(4, 4)
L, U, p = FNC.plufact(A)
A[p, :] - L @ U   # should be ≈ 0
```

Given a vector $\mathbf{b}$, we solve $\mathbf{A}\mathbf{x}=\mathbf{b}$ by first permuting the entries of $\mathbf{b}$ and then proceeding as before.

```{code-cell}
b = random.randn(4)
z = FNC.forwardsub(L, b[p])
x = FNC.backsub(U, z)
```

A residual check is successful:

```{code-cell}
b - A @ x
```

::::

While {numref}`Function {number} <function-plufact>` is a serviceable implementation, it is not the gold standard for solving linear systems in practice. Each language offers its own built-in function for PLU factorization.

::::{prf:example} Built-in PLU factorization
:label: demo-pivoting-builtin

In `scipy.solve`, the matrix `A` is PLU-factored, followed by two triangular solves. If we want to do those steps seamlessly, we can use `lu_factor` and `lu_solve`.

```{code-cell}
from scipy.linalg import lu_factor, lu_solve
A = random.randn(500, 500) 
b = ones(500)  
LU, perm = lu_factor(A)
x = lu_solve((LU, perm), b)
print(linalg.norm(b - A @ x))
```

Why would we ever bother with this? In {numref}`section-linsys-efficiency` we showed that the factorization is by far the most costly part of the solution process. A factorization object allows us to do that costly step only once per matrix, but solve with multiple right-hand sides.

```{code-cell}
start = timer()
for k in range(50): linalg.solve(A, random.rand(500))
print(f"elapsed time for 50 full solves: {timer() - start}")

start = timer()
LU, perm = lu_factor(A)
for k in range(50): lu_solve((LU, perm), random.rand(500))
print(f"elapsed time for 50 shortcut solves: {timer() - start}")
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

We construct a linear system for this matrix with $\epsilon=10^{-12}$ and exact solution $[1,1]$:

```{code-cell}
ep = 1e-12
A = array([[-ep, 1], [1, -1]])
b = A @ array([1, 1])
```

We can factor the matrix without pivoting and solve for $\mathbf{x}$.

```{code-cell}
L, U = FNC.lufact(A)
print(FNC.backsub( U, FNC.forwardsub(L, b) ))
```

Note that we have obtained only about 5 accurate digits for $x_1$. We could make the result even more inaccurate by making $\epsilon$ even smaller:

```{code-cell}
ep = 1e-20;
A = array([[-ep, 1], [1, -1]])
b = A @ array([1, 1])
L, U = FNC.lufact(A)
print(FNC.backsub( U, FNC.forwardsub(L, b) ))
```

This effect is not due to ill conditioning of the problem—a solution with PLU factorization works perfectly:

```{code-cell}
from scipy import linalg
print(linalg.solve(A, b))
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

``` python
from scipy.linalg import lu, solve
p, L, U = lu(A, p_indices=True)
x = solve(U, solve(L, b[p]))
```

Suppose instead you replace the last line above with

``` python
x = solve(solve(U, L), b[p])
```

Mathematically in terms of $\mathbf{L}$, $\mathbf{U}$, $\mathbf{p}$, and $\mathbf{b}$, what vector is found?
``````

``````{exercise}
:label: problem-pivoting-flip

✍ Suppose that `A` is a $4\times 6$ matrix and you define


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
