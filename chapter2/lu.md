---
numbering:
  enumerator: 2.4.%s
---
(section-linsys-lu)=
# LU factorization

A major tool in numerical linear algebra is to factor a given matrix into terms that are individually easier to deal with than the original. In this section we derive a means to express a square matrix using triangular factors, which will allow us to solve a linear system using forward and backward substitution.

## Outer products

```{index} outer product
```

Our derivation of the factorization hinges on an expression of matrix products in terms of vector outer products. If  $\mathbf{u}\in\real^m$ and $\mathbf{v}\in\real^n$, then the **outer product** of these vectors is the $m\times n$ matrix

:::{math}
:label: outerproduct
\mathbf{u} \mathbf{v}^T =
\begin{bmatrix}
u_1 v_1 & u_1 v_2 & \cdots & u_1 v_n \\u_2 v_1 & u_2 v_2 & \cdots & u_2 v_n \\ \vdots & \vdots & & \vdots \\ u_m v_1 & u_m v_2 & \cdots & u_m v_n
\end{bmatrix}.
:::

We illustrate the connection of outer products to matrix multiplication by a small example.

(example-lu-outer)=
::::{prf:example}
According to the usual definition of matrix multiplication,

\begin{align*}
	\small
	\begin{bmatrix}
		4 & -1  \\ -3 & 5 \\ -2 &  6	  
	\end{bmatrix}
	\begin{bmatrix}
		2 & -7 \\ -3 & 5 	  
	\end{bmatrix}
	& = 
	\small
	\begin{bmatrix}
		(4)(2) + (-1)(-3)  &  (4)(-7) + (-1)(5)   \\ 
		(-3)(2) + (5)(-3)  &  (-3)(-7) + (5)(5)  \\ 
		(-2)(2) + (6)(-3)  &  (-2)(-7) + (6)(5)  
	\end{bmatrix}.  
\end{align*}

If we break this up into the sum of two matrices, however, each is an outer product.

\begin{align*}
	& = 
	\small
	\begin{bmatrix}
		(4)(2)   &  (4)(-7)    \\ 
		(-3)(2)   &  (-3)(-7)    \\ 
		(-2)(2) &  (-2)(-7) 
	\end{bmatrix} + 
	\begin{bmatrix}
		(-1)(-3)  &  (-1)(5)  \\ 
		(5)(-3)  &  (5)(5)  \\ 
		(6)(-3)  &  (6)(5) 
	\end{bmatrix}\\[2mm]
	& = 
	\small
	\begin{bmatrix}
		4 \\ -3 \\ -2 
	\end{bmatrix} 
	\begin{bmatrix}
		2 & -7 
	\end{bmatrix} \: + \:
	\begin{bmatrix}
		-1 \\ 5 \\ 6 
	\end{bmatrix} 
	\begin{bmatrix}
		-3 & 5 
	\end{bmatrix}.
\end{align*}

Note that the vectors here are columns of the left-hand matrix and rows of the right-hand matrix. The matrix product is defined only if there are equal numbers of these.
::::

It is not hard to derive the following generalization of {numref}`Example {number} <example-lu-outer>` to all matrix products. 

::::{prf:theorem} Matrix multiplication by outer products
Write the columns of $\mathbf{A}$ as $\mathbf{a}_1,\dots,\mathbf{a}_n$ and the rows of $\mathbf{B}$ as $\mathbf{b}_1^T,\dots,\mathbf{b}_n^T$. Then

```{math}
:label: matrixouter
\mathbf{A}\mathbf{B} = \sum_{k=1}^n \mathbf{a}_k \mathbf{b}_k^T.
```
::::

## Triangular product

Equation {eq}`matrixouter` has some interesting structure for the product $\mathbf{L}\mathbf{U}$, where $\mathbf{L}$ is $n\times n$ and **lower triangular** (i.e., zero above the main diagonal) and $\mathbf{U}$ is $n\times n$ and **upper triangular** (zero below the diagonal). 

(demo-lu-outertri)=
::::{prf:example}
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #demo-lu-outertri-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-lu-outertri-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-lu-outertri-python
:::
```` 
`````
::::

Let the columns of $\mathbf{L}$ be written as $\boldsymbol{\ell}_k$ and the rows of $\mathbf{U}$ be written as $\mathbf{u}_k^T$. Then the first row of $\mathbf{L}\mathbf{U}$ is
 
```{math}
:label: outer-row1
\mathbf{e}_1^T \sum_{k=1}^n  \boldsymbol{ℓ}_k \mathbf{u}_k^T = \sum_{k=1}^n (\mathbf{e}_1^T \boldsymbol{\ell}_k) \mathbf{u}_k^T = L_{11} \mathbf{u}_1^T.
```

Likewise, the first column of $\mathbf{L}\mathbf{U}$ is

```{math}
:label: outer-col1
\left( \sum_{k=1}^n \mathbf{ℓ}_k \mathbf{u}_k^T\right) \mathbf{e}_1 = \sum_{k=1}^n \mathbf{\ell}_k (\mathbf{u}_k^T \mathbf{e}_1) = U_{11}\boldsymbol{\ell}_1.
```

These two calculations are enough to derive one of the most important algorithms in scientific computing.

## Triangular factorization

```{index} ! unit lower triangular matrix
```

Our goal is to factor a given $n\times n$ matrix $\mathbf{A}$ as the triangular product $\mathbf{A}=\mathbf{L}\mathbf{U}$. It turns out that we have $n^2+n$ total nonzero unknowns in the two triangular matrices, so we set $L_{11}=\cdots = L_{nn}=1$, making $\mathbf{L}$ a **unit lower triangular** matrix.

(demo-lu-derive)=
::::{prf:example}
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #demo-lu-derive-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-lu-derive-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-lu-derive-python
:::
```` 
`````
::::

We have arrived at the linchpin of solving linear systems. 

```{index} ! matrix factorization; LU, ! LU factorization
```
 
::::{prf:definition} LU factorization
Given $n\times n$ matrix $\mathbf{A}$, its **LU factorization** is

:::{math}
:label: def-lu
\mathbf{A} = \mathbf{L}\mathbf{U},
:::

where $\mathbf{L}$ is a unit lower triangular matrix and $\mathbf{U}$ is an upper triangular matrix.
::::

The outer product algorithm for LU factorization seen in {numref}`Demo {number} <demo-lu-derive>` is coded as {numref}`Function {number} <function-lufact>`. 

(function-lufact)=
``````{prf:algorithm} lufact
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #function-lufact-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-lufact-julia
:::
```` 


````{tab-item} Python
:sync: python
:::{embed} #function-lufact-python
:::
````
`````
``````


## Gaussian elimination and linear systems

```{index} Gaussian elimination
```

In your first matrix algebra course, you probably learned a triangularization technique called **Gaussian elimination** or row elimination to solve a linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$. In most presentations, you form an augmented matrix $[\mathbf{A}\;\mathbf{b}]$ and do row operations until the system reaches an upper triangular form, followed by backward substitution. LU factorization is equivalent to Gaussian elimination in which no row swaps are performed, and the elimination procedure produces the factors if you keep track of the row multipliers appropriately. 

Like Gaussian elimination, the primary use of LU factorization is to solve a linear system. It reduces a given linear system to two triangular ones. From this, solving $\mathbf{A}\mathbf{x}=\mathbf{b}$ follows immediately from associativity:

$$
\mathbf{b} = \mathbf{A} \mathbf{x} = (\mathbf{L} \mathbf{U}) \mathbf{x} = \mathbf{L} (\mathbf{U} \mathbf{x}).
$$

Defining $\mathbf{z} = \mathbf{U} \mathbf{x}$ leads to the following.

(algorithm-lu-solve)=
::::{prf:algorithm} Solution of linear systems by LU factorization (unstable)
1. Factor $\mathbf{L}\mathbf{U}=\mathbf{A}$.
2. Solve $\mathbf{L}\mathbf{z}=\mathbf{b}$ for $\mathbf{z}$ using forward substitution.
3. Solve $\mathbf{U}\mathbf{x}=\mathbf{z}$ for $\mathbf{x}$ using backward substitution.
::::

A key advantage of the factorization point of view is that it depends only on the matrix $\mathbf{A}$. If systems are to be solved for a single $\mathbf{A}$ but multiple different versions of $\mathbf{b}$, then the factorization approach is more efficient, as we'll see in {numref}`section-linsys-efficiency`. 

(demo-lu-solve)=
::::{prf:example}
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #demo-lu-solve-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-lu-solve-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-lu-solve-python
:::
```` 
`````
::::

As noted in the descriptions of {numref}`Function {number} <function-lufact>` and {numref}`Algorithm {number} <algorithm-lu-solve>`, the LU factorization as we have seen it so far is not stable for all matrices. In fact, it does not always even exist. The missing element is the row swapping allowed in Gaussian elimination. We will address these issues in {numref}`section-linsys-pivoting`.

## Exercises

1. ✍ For each matrix, produce an LU factorization by hand. 

    **(a)** $\quad \displaystyle \begin{bmatrix}
    2 & 3 & 4 \\
    4 & 5 & 10 \\
    4 & 8 & 2
    \end{bmatrix}\qquad$
    **(b)** $\quad \displaystyle \begin{bmatrix}
    6 & -2 & -4 & 4\\
    3 & -3 & -6 & 1 \\
    -12 & 8 & 21 & -8 \\
    -6 & 0 & -10 & 7
    \end{bmatrix}$

2. ⌨ The matrices
  
    ```{math}
    \mathbf{T}(x,y) = \begin{bmatrix}
      1 & 0 & 0 \\ 0 & 1 & 0 \\ x & y & 1
    \end{bmatrix},\qquad
    \mathbf{R}(\theta) = \begin{bmatrix}
      \cos\theta & \sin \theta & 0 \\ -\sin\theta & \cos \theta & 0 \\ 0 & 0 & 1
    \end{bmatrix}
    ```

    are used to represent translations and rotations of plane points in computer graphics. For the following, let
  
    ```{math}
    \mathbf{A} = \mathbf{T}(3,-1)\mathbf{R}(\pi/5)\mathbf{T}(-3,1), \qquad \mathbf{z} = \begin{bmatrix}
      2 \\ 2 \\ 1
    \end{bmatrix}.
    ```

    **(a)** Find $\mathbf{b} = \mathbf{A}\mathbf{z}$.

    **(b)** Use {numref}`Function {number} <function-lufact>` to find the LU factorization of $\mathbf{A}$.

    **(c)** Use the factors with triangular substitutions in order to solve $\mathbf{A}\mathbf{x}=\mathbf{b}$, and find $\mathbf{x}-\mathbf{z}$.
  
    (problem-bigcorner)=
3. ⌨ Define
  
    ```{math}
    \mathbf{A}= \begin{bmatrix}
      1 & 0 & 0 & 0 & 10^{12} \\
      1 & 1 & 0 & 0 & 0 \\
      0 & 1 & 1 & 0 & 0 \\
      0 & 0 & 1 & 1 & 0 \\
      0 & 0 & 0 & 1 & 0
    \end{bmatrix},
    \quad \hat{\mathbf{x}} = \begin{bmatrix}
      0 \\ 1/3 \\ 2/3 \\ 1 \\ 4/3
    \end{bmatrix},
    \quad \mathbf{b} = \mathbf{A}\hat{\mathbf{x}}.
    ```

    **(a)** Using {numref}`Function {number}<function-lufact>` and triangular substitutions, solve the linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$, showing the result. To the nearest integer, how many accurate digits are in the result? (The answer is much less than the full 16 of double precision.)

    **(b)** Repeat part (a) with $10^{20}$ as the element in the upper right corner. (The result is even less accurate. We will study the causes of such low accuracy in {numref}`section-linsys-condition-number`.)
  
4. ⌨ Let

    $$
	\mathbf{A} = 
	\begin{bmatrix}
     1 & 1 & 0 & 1 & 0 & 0 \\
     0 & 1 & 1 & 0 & 1 & 0 \\
     0 & 0 & 1 & 1 & 0 & 1 \\
     1 & 0 & 0 & 1 & 1 & 0 \\
     1 & 1 & 0 & 0 & 1 & 1 \\
     0 & 1 & 1 & 0 & 0 & 1
    \end{bmatrix}.
    $$

    Verify computationally that if $\mathbf{A}=\mathbf{L}\mathbf{U}$ is the LU factorization, then the elements of $\mathbf{L}$, $\mathbf{U}$, $\mathbf{L}^{-1}$, and $\mathbf{U}^{-1}$ are all integers. Do **not** rely just on visual inspection of the numbers; perform a more definitive test.

5. ⌨ {numref}`Function {number}<function-lufact>` factors $\mathbf{A}=\mathbf{L}\mathbf{U}$ in such a way that $\mathbf{L}$ is a unit lower triangular matrix—that is, has all ones on the diagonal. It is also possible to define the factorization so that $\mathbf{U}$ is a unit upper triangular matrix instead. Write a function `lufact2` that uses {numref}`Function {number}<function-lufact>` *without modification* to produce this version of the factorization. (Hint: Begin with the standard LU factorization of $\mathbf{A}^T$.) Demonstrate on a nontrivial $4\times 4$ example.

6. When computing the determinant of a matrix by hand, it's common to use cofactor expansion and apply the definition recursively. But this is terribly inefficient as a function of the matrix size.
  
    **(a)** ✍ Explain using determinant properties why, if $\mathbf{A}=\mathbf{L}\mathbf{U}$ is an LU factorization,

    ```{math}
      \det(\mathbf{A}) = U_{11}U_{22}\cdots U_{nn}=\prod_{i=1}^n U_{ii}.
    ```

    **(b)** ⌨ Using the result of part (a), write a function `determinant(A)` that computes the determinant using {numref}`Function {number}<function-lufact>`. Test your function on at least two nontriangular $5\times 5$ matrices, comparing your result to the result of the standard `det` function.
