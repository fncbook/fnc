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

## Vectors and matrices in Julia

In Julia, vectors and matrices are one-dimensional and two-dimensional arrays, respectively. A lot of how Julia deals with them is easy to remember once learned. There's a lot to learn, though, so we give just some of the basics here, and we will pick up more as we go from the code used in our examples and functions.  

(demo-matrices-julia)=
```{prf:example}
```





:::{index} ! Julia; size, ! Julia; length
:::

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item}
:columns: 7


Square brackets are used to enclose elements of a matrix or vector. Use spaces for horizontal concatenation, and semicolons or new lines to indicate vertical concatenation.


:::
:::{grid-item-card}
:columns: 5



The `size` function returns the number of rows and columns in a matrix. Use `length` to get the number of elements in a vector or matrix.

:::
::::


```{code-cell}
A = [ 1 2 3 4 5; 50 40 30 20 10
    π sqrt(2) exp(1) (1+sqrt(5))/2 log(3) ]
```

```{code-cell}
m,n = size(A)
```
::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item}
:columns: 7


A vector is not quite the same thing as a matrix. It has only one dimension, not two. Separate its elements by commas or semicolons.


:::
:::{grid-item-card}
:columns: 5


Separate elements of a vector by commas or semicolons. If you use spaces, you will get a matrix with one row, which is treated differently.

:::
::::

```{code-cell}
x = [ 3, 3, 0, 1, 0 ]
size(x)
```

For many purposes, however, an $n$-vector in Julia is a lot like an $n\times 1$ column vector.

Concatenated elements within brackets may be matrices or vectors for a block representation, as long as all the block sizes are compatible.

```{code-cell}
[ x  x ]
```

```{code-cell}
[ x; x ]
```

The `zeros` and `ones` functions construct matrices with entries all zero or one, respectively.

```{code-cell}
B = [ zeros(3,2) ones(3,1) ]
```

```{index} ! Julia; transpose, ! Julia; adjoint, ! Julia; \'
```

A single quote `'` after a matrix returns its adjoint. For real matrices, this is the transpose; for complex-valued matrices, the elements are also conjugated. 

```{code-cell}
A'
```

If `x` is simply a vector, then its transpose has a row shape.

```{code-cell}
x'
```

```{index} ! Julia; range, ! Julia; \:
```

There are many convenient shorthand ways of building vectors and matrices other than entering all of their entries directly or in a loop. To get a range with evenly spaced entries between two endpoints, you have two options. One is to use a colon `:`.

```{code-cell}
y = 1:4              # start:stop
```

```{code-cell}
z = 0:3:12           # start:step:stop
```

(Ranges are not strictly considered vectors, but they behave identically in most circumstances.) Instead of specifying the step size, you can give the number of points in the range if you use `range`.

```{code-cell}
s = range(-1,1,length=5)
```

:::{index} ! Julia; end, ! Julia; indexing arrays
:::

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item}
:columns: 7


Accessing an element is done by giving one (for a vector) or two (for a matrix) index values within square brackets. 


:::
:::{grid-item-card}
:columns: 5


The `end` keyword refers to the last element in a dimension. It saves you from having to compute and store the size of the matrix first.

:::
::::

```{code-cell}
a = A[2,end-1]
```

```{code-cell}
x[2]
```

The indices can be vectors or ranges, in which case a block of the matrix is accessed.

```{code-cell}
A[1:2,end-2:end]    # first two rows, last three columns
```

```{index} Julia; \:
```

If a dimension has only the index `:` (a colon), then it refers to all the entries in that dimension of the matrix.

```{code-cell}
A[:,1:2:end]        # all of the odd columns
```

:::{index} ! Julia; diagm
:::

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item}
:columns: 7


The matrix and vector senses of addition, subtraction, scalar multiplication, multiplication, and power are all handled by the usual symbols. 


:::
:::{grid-item-card}
:columns: 5


Use `diagm` to construct a matrix by its diagonals. A more general syntax puts elements on super- or subdiagonals.

:::
::::

```{code-cell}
B = diagm( [-1,0,-5] )   # create a diagonal matrix
```

```{code-cell}
BA = B*A     # matrix product
```
::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item}
:columns: 7


`A*B` causes an error because the dimensions aren't compatible.


:::
:::{grid-item-card}
:columns: 5


Errors are formally called *exceptions* in Julia.

:::
::::

```{code-cell} julia
A*B    # throws an error
```

A square matrix raised to an integer power is the same as repeated matrix multiplication.

```{code-cell}
B^3    # same as B*B*B
```

Sometimes one instead wants to treat a matrix or vector as a mere array and simply apply a single operation to each element of it. For multiplication, division, and power, the corresponding operators start with a dot.

```{code-cell}
C = -A;
```

`A*C` would be an error.

```{code-cell}
elementwise = A.*C
```

```{index} Julia; broadcasting
```

The two operands of a dot operator have to have the same size—unless one is a scalar, in which case it is expanded or *broadcast* to be the same size as the other operand.

```{code-cell}
x_to_two = x.^2
```

```{code-cell}
two_to_x = 2 .^ x
```

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item}
:columns: 7


Most of the mathematical functions, such as cos, sin, log, exp, and sqrt, expect scalars as operands. However, you can broadcast any function, including ones that you have defined, across a vector or array by using a special dot syntax.


:::
:::{grid-item-card}
:columns: 5


A dot added to the end of a function name means to apply the function elementwise to an array.

:::
::::

```{code-cell}
show(cos.(π*x))    # broadcast to a function
```

Rather than dotting multiple individual functions, you can use `@.` before an expression to broadcast everything within it.

```{code-cell}
show(@. cos(π*(x+1)^3))    # broadcast an entire expression
```





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

Being able to extract specific rows and columns of a matrix makes it straightforward to do row- and column-oriented operations, such as linear combinations.

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

```{index} ! Julia; +=
```

The Julia equivalent is

```julia
A[:,1] += 2A[:,3]
```

The `+=` operator means to increment the item on the left-hand side. There are similar interpretations for `-=` and `*=`.
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
