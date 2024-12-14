---
numbering:
  enumerator: 2.2.%s
---
(section-linsys-matrices)=
# Computing with matrices

```{attention}
We recommend that you review the linear algebra material in {doc}`../appendix/linear-algebra` before reading this section.
```

At a reductive level, a matrix is a table of numbers that obeys certain algebraic laws. But matrices are pervasive in scientific computation, mainly because they represent linear operations on vectors. Moreover, vectors go far beyond the three-dimensional representations of physical quantities you learned about in calculus. 

## Notation

We use capital letters in bold to refer to matrices, and lowercase bold letters for vectors. *All named vectors in this book are column vectors*. The bold symbol $\boldsymbol{0}$ may refer to a vector of all zeros or to a zero matrix, depending on context; we use $0$ as the scalar zero only.

To refer to a specific element of a matrix, we use the uppercase name of the matrix *without* boldface, as in $A_{24}$ to mean the $(2,4)$ element of $\mathbf{A}$.[^subscript] To refer to an element of a vector, we use just one subscript, as in $x_3$. If you see a boldface character with one or more subscripts, then you know that it is a matrix or vector that belongs to a sequence or indexed collection.

[^subscript]: This aspect of our notation is slightly unusual. More frequently one would see the lowercase $a_{24}$ in this context. We feel that our notation lends more consistency and clarity to expressions with mixed symbols, and it is more like how computer code is written.

We will have frequent need to refer to the individual columns of a matrix as vectors. Our convention is to use a lowercase bold version of the matrix name with a subscript to represent the column number. Thus, $\mathbf{a}_1,\mathbf{a}_2,\ldots,\mathbf{a}_n$ are the columns of the $m\times n$ matrix $\mathbf{A}$. Conversely, whenever we define a sequence of vectors $\mathbf{v}_1,\ldots,\mathbf{v}_p$, we can implicitly consider them to be columns of a matrix $\mathbf{V}$. Sometimes we might write $\mathbf{V}=\bigl[ \mathbf{v}_j \bigr]$ to emphasize the connection.

```{index} transpose, adjoint of a matrix, ! symmetric matrix, hermitian matrix
```

The notation $\mathbf{A}^T$ is used for the transpose of a matrix, whether it is real or complex. In the case of complex matrices, it's almost always more desirable to use the **adjoint** $\mathbf{A}^*$, which is the transpose with the complex conjugate of each element.[^conjugate]  If $\mathbf{A}$ is real, then $\mathbf{A}^*=\mathbf{A}^T$. A **symmetric matrix** is a square matrix such that $\mathbf{A}^T=\mathbf{A}$. 

[^conjugate]: The conjugate of a complex number is found by replacing all references to the imaginary unit $i$ by $-i$.

```{index} identity matrix
```

The **identity matrix** of size $n$ is denoted $\mathbf{I}$, or sometimes $\mathbf{I}_n$ if emphasizing the size is important in context. For columns of the identity we break with our usual naming convention and denote them by $\mathbf{e}_j$.

## Block matrix expressions

We will often find it useful to break a matrix into separately named pieces. For example, we might write

```{math}
  \mathbf{A} =
  \begin{bmatrix}
    \mathbf{A}_{11} & \mathbf{A}_{12} & \mathbf{A}_{13} \\
    \mathbf{A}_{21} & \mathbf{A}_{22} & \mathbf{A}_{23}
  \end{bmatrix}, \qquad
  \mathbf{B} =
  \begin{bmatrix}
    \mathbf{B}_1 \\ \mathbf{B}_2 \\ \mathbf{B}_3
  \end{bmatrix}.
```

It's understood that blocks that are on top of one another have the same number of columns, and blocks that are side by side have the same number of rows. Typically, if the blocks all have compatible dimensions, then they can be multiplied as though the blocks were scalars. For instance, continuing with the definitions above, we say that $\mathbf{A}$ is block-$2\times 3$ and $\mathbf{B}$ is block-$3\times 1$, so we can write

```{math}
  \mathbf{A} \mathbf{B} =
  \begin{bmatrix}
    \mathbf{A}_{11}\mathbf{B}_1 + \mathbf{A}_{12}\mathbf{B}_2 + \mathbf{A}_{13}\mathbf{B}_3 \\
    \mathbf{A}_{21}\mathbf{B}_1 + \mathbf{A}_{22}\mathbf{B}_2 + \mathbf{A}_{23}\mathbf{B}_3
  \end{bmatrix},
```

provided that the individual block products are well-defined. For transposes we have, for example,

```{math}
  \mathbf{A}^T =
  \begin{bmatrix}
    \mathbf{A}_{11}^T & \mathbf{A}_{21}^T \\[2mm]
    \mathbf{A}_{12}^T & \mathbf{A}_{22}^T \\[2mm]
    \mathbf{A}_{13}^T & \mathbf{A}_{23}^T
  \end{bmatrix}.
```

## Vector and matrix basics 

Vectors and matrices are integral to scientific computing. All modern languages provide ways to work with them beyond manipulation of individual elements. 

(demo-matrices)=
::::{prf:example}
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #demo-matrices-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-matrices-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-matrices-python
:::
```` 
`````
::::



## Row and column operations

A critical identity in matrix multiplication is

```{math}
  \mathbf{A} \mathbf{e}_j = \mathbf{a}_j.
```

```{prf:observation}
Multiplication on the right by column $j$ of the identity reproduces the $j$th column of a matrix. 
```

Furthermore, the expression

```{math}
  \mathbf{A}
  \begin{bmatrix}
    \mathbf{e}_1 & \mathbf{e}_3 & \mathbf{e}_5
  \end{bmatrix}
```

reproduces three columns. An equivalent expression in Julia would be `A[:,1:2:5]`. 

We can extend the same idea to rows by using the general identity $(\mathbf{R}\mathbf{S})^T=\mathbf{S}^T\mathbf{R}^T$. Let $\mathbf{B}=\mathbf{A}^T$ have columns $\bigl[ \mathbf{b}_j \bigr]$, and note

```{math}
  (\mathbf{b}_j)^T = (\mathbf{B} \mathbf{e}_j)^T = \mathbf{e}_j^T \mathbf{B}^T = \mathbf{e}_j^T \mathbf{A}.
```

But $\mathbf{e}_j^T$ is the $j$th *row* of $\mathbf{I}$, and $\mathbf{b}_j^T$ is the transpose of the $j$th column of $\mathbf{B}$, which is the $j$th *row* of $\mathbf{A}$ by $\mathbf{B}=\mathbf{A}^T$. Thus, *multiplication on the left by row $j$ of the identity extracts the $j$th row*. Extracting the single element $(i,j)$ from the matrix is, therefore, $\mathbf{e}_i^T \mathbf{A} \mathbf{e}_j$.

Being able to extract specific rows and columns of a matrix via algebra makes it straightforward to do row- and column-oriented operations, such as linear combinations.

````{prf:example}
 Say that $\mathbf{A}$ has five columns. Adding twice the third column of $\mathbf{A}$ to its first column is done by

```{math}
\mathbf{A}(\mathbf{e}_1+2\mathbf{e}_3).
```

Suppose we want to do this operation "in place," meaning replacing the first column of $\mathbf{A}$ with this value and leaving the other four columns of $\mathbf{A}$ alone. We can replace $\mathbf{A}$ with

```{math}
  \mathbf{A}
  \begin{bmatrix}
    \mathbf{e}_1+2\mathbf{e}_3 & \mathbf{e}_2 & \mathbf{e}_3 & \mathbf{e}_4 & \mathbf{e}_5
  \end{bmatrix}.
```

:::::{tab-set}
::::{tab-item} Julia
:sync: julia

```{index} ! Julia; +=
```

The Julia equivalent is

```julia
A[:, 1] += 2A[:, 3]
```

The `+=` operator means to increment the item on the left-hand side. There are similar interpretations for `-=` and `*=`.
::::

:::{tab-item} MATLAB
:sync: matlab
The MATLAB equivalent is

```matlab
A(:, 1) = A(:, 1) + 2 * A(:, 3)
```

:::
:::{tab-item} Python
:sync: python
The NumPy equivalent is

```python
A[:, 0] += 2 * A[:, 2]
```

The `+=` operator means to increment the item on the left-hand side. There are similar interpretations for `-=` and `*=`
:::
::::

````

## Exercises

1. ✍ Suppose 
   
    $$
    \mathbf{C} =
    \begin{bmatrix}
      \mathbf{I} & \mathbf{A} \\ -\mathbf{I} & \mathbf{B}
    \end{bmatrix}.
    $$ 
    
    Using block notation, find $\mathbf{C}^2$ and $\mathbf{C}^3$.

2. ⌨  Let

    ```{math}
    \mathbf{A} =
    \begin{bmatrix}
      2 & 1 & 1 & 0 \\ 0 & -1 & 4 & 1 \\ 2 & 2 & 0 & -2 \\ 1 & 3 & -1
      & 5
    \end{bmatrix},
    \quad
    \mathbf{B} =
    \begin{bmatrix}
      3 & -1 & 0 & 2 \\ 7 & 1 & 0 & 2
    \end{bmatrix},
    ```

    ```{math}
    \mathbf{u} =
    \begin{bmatrix}
      2 \\ -1 \\ 3 \\ 1
    \end{bmatrix},
    \quad
    \mathbf{v} =
    \begin{bmatrix}
      \pi \\ e
    \end{bmatrix}.
    ```

    (Do not round off the values in $\mathbf{v}$—find them using native Julia commands.) For each expression below, use Julia to find the result, or explain why the result does not exist.

    **(a)** $\mathbf{A}\mathbf{B},\quad$
    **(b)** $\mathbf{B} \mathbf{A},\quad$
    **(c)** $\mathbf{v}^T \mathbf{B},\quad$
    **(d)** $\mathbf{B} \mathbf{u},\quad$
    **(e)** $\bigl[ \, \mathbf{u}\:\: \mathbf{A}\mathbf{u} \:\: \mathbf{A}^2 \mathbf{u} \:\: \mathbf{A}^3 \mathbf{u} \bigr]$.

3. ⌨  Let
  
    ```{math}
    \mathbf{u} =
    \begin{bmatrix}
      1\\3\\5\\7\\9\\11
    \end{bmatrix}, \qquad
    \mathbf{v} =
    \begin{bmatrix}
      -60 \\ -50 \\ -40 \\ -30 \\ -20 \\ -10
    \end{bmatrix}.
    ```

    Find the inner products $\mathbf{u}^T\mathbf{v}$ and $\mathbf{v}^T\mathbf{u}$ and the outer products $\mathbf{u}\mathbf{v}^T$ and $\mathbf{v}\mathbf{u}^T$.

4. ⌨ In Julia, give a demonstration of the identity $(\mathbf{A}\mathbf{B})^T=\mathbf{B}^T\mathbf{A}^T$ for some arbitrarily chosen $3\times 4$ matrix $\mathbf{A}$ and $4\times 2$ matrix $\mathbf{B}$.

   (problem-inverseprod)=
5. ✍ Prove that if $\mathbf{A}$ and $\mathbf{B}$ are invertible, then $(\mathbf{A}\mathbf{B})^{-1}=\mathbf{B}^{-1}\mathbf{A}^{-1}$. (In producing the inverse, it follows that $\mathbf{A}\mathbf{B}$ is invertible as well.)

6. ✍ Suppose $\mathbf{B}$ is an arbitrary $4\times 3$ matrix. In each part below a matrix $\mathbf{A}$ is described in terms of $\mathbf{B}$. Express $\mathbf{A}$ as a product of $\mathbf{B}$ with one or more other matrices.
  
    **(a)** $\mathbf{A}\in\mathbb{R}^{4 \times 1}$ is the result of adding the first column of $\mathbf{B}$ to $-2$ times the last column of $\mathbf{B}$.

    **(b)** The rows of $\mathbf{A}\in\mathbb{R}^{4 \times 3}$ are the rows of $\mathbf{B}$ in order 4,3,2,1.

    **(c)** The first column of $\mathbf{A}\in\mathbb{R}^{4 \times 3}$ is $1$ times the first column of $\mathbf{B}$, the second column of $\mathbf{A}$ is $2$ times the second column of $\mathbf{B}$,
    and the third column of $\mathbf{A}$ is $3$ times the third column of $\mathbf{B}$.

    **(d)** $A$ is the scalar sum of all elements of $\mathbf{B}$.
  
7.  **(a)** ✍ Prove that for real vectors $\mathbf{v}$ and $\mathbf{w}$ of the same length, the inner products $\mathbf{v}^T\mathbf{w}$ and $\mathbf{w}^T\mathbf{v}$ are equal.

    **(b)** ✍ Prove true, or give a counterexample for, the equivalent statement about outer products, $\mathbf{v}\mathbf{w}^T$ and $\mathbf{w}\mathbf{v}^T$.
