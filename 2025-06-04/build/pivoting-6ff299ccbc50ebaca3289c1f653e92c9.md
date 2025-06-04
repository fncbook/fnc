---
numbering:
  enumerator: 2.6.%s
---
(section-linsys-pivoting)=

# Row pivoting

```{index} matrix factorization; LU
```

As mentioned in {numref}`section-linsys-lu`, the $\mathbf{A}=\mathbf{L}\mathbf{U}$ factorization is not stable for every nonsingular $\mathbf{A}$. Indeed, the factorization does not always even exist.

::::{prf:example} Failure of naive LU factorization
:label: demo-pivoting-fail

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-pivoting-fail-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-pivoting-fail-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-pivoting-fail-python
:::
```` 
`````

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

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-pivoting-fix-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-pivoting-fix-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-pivoting-fix-python
:::
```` 
`````

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

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-pivoting-permute-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-pivoting-permute-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-pivoting-permute-python
:::
```` 
`````

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

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #function-plufact-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-plufact-julia
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #function-plufact-python
:::
````
`````
``````

Ideally, the PLU factorization takes $\sim \frac{2}{3}n^3$ flops asymptotically, just like LU without pivoting. The implementation in {numref}`Function {number} <function-plufact>` does not achieve this optimal flop count, however. Like {numref}`Function {number} <function-lufact>`, it does unnecessary operations on structurally known zeros for the sake of being easier to understand.

## Linear systems

The output of {numref}`Function {number} <function-plufact>` is a factorization of a row-permuted $\mathbf{A}$. Therefore, given a linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$, we have to permute $\mathbf{b}$ the same way before applying forward and backward substitution. This is equivalent to changing the order of the equations in a linear system, which does not affect its solution.

::::{prf:example} PLU factorization for solving linear systems
:label: demo-pivoting-usage

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-pivoting-usage-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-pivoting-usage-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-pivoting-usage-python
:::
```` 
`````

::::

The `lu` function from the built-in package `LinearAlgebra` returns the same three outputs as {numref}`Function {number} <function-plufact>`. If you only request one output, it will be a factorization object that can be used with a backslash. This is useful when you want to solve with multiple versions of $\mathbf{b}$ but do the factorization only once.

::::{prf:example} Built-in PLU factorization
:label: demo-pivoting-builtin

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-pivoting-builtin-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-pivoting-builtin-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-pivoting-builtin-python
:::
```` 
`````

::::

## Stability

```{index} stability
```

There is one detail of the row pivoting algorithm that might seem arbitrary: why choose the pivot of largest magnitude in a column, rather than, say, the uppermost nonzero in the column? The answer is numerical stability.

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

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-pivoting-stable-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-pivoting-stable-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-pivoting-stable-python
:::
```` 
`````

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

✍ Let $\mathbf{A}$ be a square matrix and $\mathbf{b}$ be a column vector of compatible length. Here is correct Julia code to solve $\mathbf{A}\mathbf{x}=\mathbf{b}$:

``` julia
L,U,p = lu(A)
x = U \ (L\b[p])
```

Suppose instead you replace the last line above with

``` julia
x = U \ L \ b[p]
```

Mathematically in terms of $\mathbf{L}$, $\mathbf{U}$, $\mathbf{p}$, and $\mathbf{b}$, what vector is found?
``````

``````{exercise}
:label: problem-pivoting-flip

✍ Suppose that `A` is a $4\times 6$ matrix in Julia and you define

``` julia
B = A[end:-1:1,end:-1:1]
```

Show that $\mathbf{B} = \mathbf{P} \mathbf{A} \mathbf{Q}$ for certain matrices $\mathbf{P}$ and $\mathbf{Q}$.
``````

``````{exercise}
:label: problem-pivoting-perminverse
✍ An $n\times n$ *permutation matrix* $\mathbf{P}$ is a reordering of the rows of an identity matrix such that $\mathbf{P} \mathbf{A}$  has the effect of moving rows $1,2,\ldots,n$ of $\mathbf{A}$ to new positions $i_1,i_2,\ldots,i_n$. Then $\mathbf{P}$ can be expressed as

```{math}
\mathbf{P} = \mathbf{e}_{i_1}\mathbf{e}_1^T + \mathbf{e}_{i_2}\mathbf{e}_2^T + \cdots + \mathbf{e}_{i_n}\mathbf{e}_n^T.
```

**(a)** For the case $n=4$ and $i_1=3$, $i_2=2$, $i_3=4$, $i_4=1$, write out separately, as matrices, all four of the terms in the sum. Then add them together to find $\mathbf{P}$.

**(b)** Use the formula in the general case to show that $\mathbf{P}^{-1}=\mathbf{P}^T$.
``````
