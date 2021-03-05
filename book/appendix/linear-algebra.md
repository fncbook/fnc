# Review: Linear algebra

```{glossary}
inverse
  Multiplicative inverse of a matrix.

invertible
  Describing a matrix (necessarily square) that has an inverse.
```

An $m\times n$ matrix $\mathbf{A}$ is a rectangular $m$-by-$n$ array of numbers called **elements** or **entries**.  The numbers $m$ and $n$ are called the **row dimension** and the **column dimension**, respectively; collectively they describe the **size** of $\mathbf{A}$. We say $\mathbf{A}$ belongs to $\mathbb{R}^{m\times n}$ if its entries are real or $\mathbb{C}^{m\times n}$ if they are complex-valued.  A **square** matrix has equal row and column dimensions. A **row vector** has dimension $1\times n$, while a **column vector** has dimension $m \times 1$. By default, a vector is understood to be a column vector, and we use $\mathbb{R}^n$ or $\mathbb{C}^n$ to denote spaces of vectors.  An ordinary number in $\mathbb{R}$ or $\mathbb{C}$ may be called a **scalar**.

We use capital letters in bold to refer to matrices, and lowercase bold letters for vectors. In this book, all vectors are column vectors—in other words, matrices with multiple rows and one column. The bold symbol $\boldsymbol{0}$ may refer to a vector of all zeros or to a zero matrix, depending on context; we use $0$ as the scalar zero only.

To refer to a specific element of a matrix, we use the uppercase name of the matrix *without* boldface, as in $A_{24}$ to mean the $(2,4)$ element of $\mathbf{A}$. To refer to an element of a vector, we use just one subscript, as in $x_3$. If you see a boldface character with one or more subscripts, then you know that it is a matrix or vector that belongs to a numbered collection.

We will have frequent need to refer to the individual columns of a matrix as vectors. Our convention, is to use a lowercase bold version of the matrix name, with a subscript to represent the column number. Thus, $\mathbf{a}_1,\mathbf{a}_2,\ldots,\mathbf{a}_n$ are the columns of the $m\times n$ matrix $\mathbf{A}$. Conversely, whenever we define a sequence of vectors $\mathbf{v}_1,\ldots,\mathbf{v}_p$, we can implicitly consider them to be columns of a matrix $\mathbf{V}$. Sometimes we might write $\mathbf{V}=\bigl[ \mathbf{v}_j \bigr]$ to emphasize the connection.

The **diagonal** (main diagonal) of an $n\times n$ matrix $\mathbf{A}$ refers to the entries $A_{ii}$, $i=1,\ldots,n$. The entries $A_{ij}$ where $j-i=k$ are on a **superdiagonal** if $k>0$ and a **subdiagonal** if $k<0$. The diagonals are numbered as suggested here:

```{math}
  \begin{bmatrix}
    0 & 1 & 2 & \cdots & n-1 \\
    -1 & 0 & 1 & \cdots & n-2 \\
    \vdots & \ddots & \ddots & \ddots & \vdots \\
    -n+2 & \cdots & -1 & 0 & 1\\
    -n+1 & \cdots & -2 & -1 & 0
  \end{bmatrix}.
```

A **diagonal** matrix is one whose entries are all zero off the main diagonal.  An **upper triangular** matrix $\mathbf{U}$ has entries $U_{ij}$ with $U_{ij}=0$ if $i>j$, and a **lower triangular** matrix $\mathbf{L}$ has $L_{ij}=0$ if $i<j$.

The **transpose** of $\mathbf{A}\in\mathbb{C}^{m\times n}$ is the matrix $\mathbf{A}^T\in\mathbb{C}^{n\times m}$ given by

```{math}
  \mathbf{A}^T =
  \begin{bmatrix}
    A_{11} & A_{21} & \cdots & A_{m1}\\
    \vdots & \vdots & & \vdots\\
    A_{1n} & A_{2n} & \cdots & A_{mn}
  \end{bmatrix}.
```

The **conjugate transpose** or **hermitian** is given by $\mathbf{A}^*=\overline{\mathbf{A}^T}$, where the bar denotes taking a complex conjugate.[^conj] If $\mathbf{A}$ is real, then $\mathbf{A}^*=\mathbf{A}^T$. A square matrix is **symmetric** if $\mathbf{A}^T=\mathbf{A}$ and **hermitian** if $\mathbf{A}^*=\mathbf{A}$.

[^conj]: The conjugate of a complex number is found by   replacing all references to the imaginary unit $i$ by $-i$. We will not see much of complex numbers until the second half of the book.

## Algebra

Matrices of the same size may be added componentwise.  Multiplication by a scalar is also defined elementwise. These operations obey the familiar laws of commutativity, associativity, and distributivity. The multiplication of two matrices, on the other hand, is much less straightforward.

In order for matrices $\mathbf{A}$ and $\mathbf{B}$ to be multiplied, it is necessary that their "inner" dimensions match—i.~e., $\mathbf{A}$ is $m\times p$ and $\mathbf{B}$ is $p \times n$. Note that even if $\mathbf{A}\mathbf{B}$ is defined, $\mathbf{B}\mathbf{A}$ may not be, unless $m=n$. In terms of scalar components, the $(i,j)$ entry of $\mathbf{C}=\mathbf{A}\mathbf{B}$ is given by

```{math}
  :label: scalarmatrixmult
  C_{ij} = \sum_{k=1}^p A_{ik} B_{kj}.
```

An important identity is that when $\mathbf{A}\mathbf{B}$ is defined,

```{math}
  :label: transposeidentity
  (\mathbf{A}\mathbf{B})^T=\mathbf{B}^T\mathbf{A}^T
```

The latter product is always defined in this situation.

```{margin}
Matrix multiplication is not commutative.
```

However, $\mathbf{A}\mathbf{B}= \mathbf{B}\mathbf{A}$ is *not* true in general, even when both products are defined. That is, matrix multiplication is not commutative. It is, however, associative: $\mathbf{A}\mathbf{B}\mathbf{C}=(\mathbf{A}\mathbf{B})\mathbf{C}=\mathbf{A}(\mathbf{B}\mathbf{C})$. Changing the ordering of operations by associativity is a trick we will use repeatedly.

It is worth reinterpreting {eq}`scalarmatrixmult` at a vector level. If $\mathbf{A}$ has dimensions $m\times n$, it can be multiplied on the right by an $n \times 1$ column vector $\mathbf{v}$ to produce an $m \times 1$ column vector $\mathbf{A}\mathbf{v}$, which satisfies

```{math}
  :label: mvcol
  \mathbf{A}\mathbf{v} =
  \begin{bmatrix}
    \displaystyle \sum_k A_{1k}v_k \\[2mm]
    \displaystyle \sum_k A_{2k}v_k \\
    \vdots\\
    \displaystyle \sum_k A_{mk}v_k
  \end{bmatrix}
  = v_1
  \begin{bmatrix}
    A_{11}\\A_{21}\\\vdots\\A_{m1}
  \end{bmatrix} +
  v_2
  \begin{bmatrix}
    A_{12}\\A_{22}\\\vdots\\A_{m2}
  \end{bmatrix} +
  \cdots + v_n
  \begin{bmatrix}
    A_{1n}\\A_{2n}\\\vdots\\A_{mn}
  \end{bmatrix} = v_1 \mathbf{a}_1 + \cdots + v_n \mathbf{a}_n.
```

In words, we say that $\mathbf{A}\mathbf{v}$ is a  **linear combination** of the columns of $\mathbf{A}$. Equation {eq}`mvcol` is very important: \emph{Multiplying a matrix on the right by a column vector produces a linear combination of the columns of the matrix}. There is a similar interpretation of multiplying $\mathbf{A}$ on the left by a row vector. Keeping to our convention that boldface letters represent *column* vectors, we write, for $\mathbf{v}\in\mathbb{R}^m$,

```{math}
  :label: mvrow
  \mathbf{v}^T \mathbf{A} =
  \begin{bmatrix}
    \displaystyle \sum_k v_k A_{k1} & \displaystyle \sum_k v_k A_{k2} & \cdots & \displaystyle \sum_k v_k A_{kn}
  \end{bmatrix}
   = v_1
  \begin{bmatrix}
    A_{11} &  \cdots & A_{1n}
  \end{bmatrix} +
  v_2
  \begin{bmatrix}
    A_{21} & \cdots & A_{2n}
  \end{bmatrix} +
  \cdots + v_m
  \begin{bmatrix}
    A_{m1} & \cdots & A_{mn}
  \end{bmatrix}.
```

Thus multiplying a matrix on the left by a row vector produces a linear combination of the rows of the matrix.

These two observations extend to more general matrix-matrix multiplications. One can show that (here again $\mathbf{A}$ is $m\times p$ and $\mathbf{B}$ is $p\times n$)

```{math}
  :label: mmhoriz
  \mathbf{A}\mathbf{B} =
  \mathbf{A} \begin{bmatrix}
    \mathbf{b}_1 & \mathbf{b}_2 & \cdots & \mathbf{b}_n
  \end{bmatrix}
  = \begin{bmatrix}
    \mathbf{A}\mathbf{b}_1 & \mathbf{A}\mathbf{b}_2 & \cdots & A\mathbf{b}_n
  \end{bmatrix}
```

In words, a matrix-matrix product is a horizontal concatenation of matrix-vector products involving the columns of the right-hand matrix. Equivalently, if we write $\mathbf{A}$ in terms of rows, then

```{math}
  :label: mmvert
  \mathbf{A} =
  \begin{bmatrix}
    \mathbf{w}_1^T \\[2mm] \mathbf{w}_2^T \\ \vdots \\ \mathbf{w}_m^T
  \end{bmatrix}
  \qquad\text{ implies }\qquad
    \mathbf{A}\mathbf{B} =
  \begin{bmatrix}
    \mathbf{w}_1^T \mathbf{B} \\[2mm] \mathbf{w}_2^T \mathbf{B} \\ \vdots \\ \mathbf{w}_m^T \mathbf{B}
  \end{bmatrix}.
```

Thus, a matrix-matrix product is also a vertical concatenation of vector-matrix products involving the rows of the left-hand matrix. All of our representations of matrix multiplication are equivalent, so whichever one is most convenient at any moment can be used.

````{prf:example}
Let
```{math}
\mathbf{A} = \begin{bmatrix}
      1 & -1 \\ 0 & 2 \\ -3 & 1
    \end{bmatrix}, \qquad
\mathbf{B} = \begin{bmatrix}
      2 & -1 & 0 & 4 \\ 1 & 1 & 3 & 2
    \end{bmatrix}.
```

Then, going by {eq}`scalarmatrixmult`, we get

```{math}
\begin{split}
  \mathbf{A}\mathbf{B} &= \begin{bmatrix}
          (1)(2) + (-1)(1) & (1)(-1) + (-1)(1) & (1)(0) + (-1)(3) & (1)(4) + (-1)(2) \\
          (0)(2) + (2)(1) & (0)(-1) + (2)(1) & (0)(0) + (2)(3) & (0)(4) + (2)(2) \\
          (-3)(2) + (1)(1) & (-3)(-1) + (1)(1) & (-3)(0) + (1)(3) & (-3)(4) + (1)(2)
        \end{bmatrix} \\
  &= \begin{bmatrix}
      1 & -2 & -3 & 2 \\ 2 & 2 & 6 & 4 \\ -5 & 4 & 3 & -10
    \end{bmatrix}.
\end{split}
```

But note also, for instance, that

```{math}
  \mathbf{A} \begin{bmatrix} 2 \\ 1 \end{bmatrix} = 2 \begin{bmatrix} 1 \\ 0 \\ -3
    \end{bmatrix} + 1 \begin{bmatrix} -1 \\ 2 \\ 1 \end{bmatrix} = \begin{bmatrix} 1 \\ 2 \\ -5 \end{bmatrix},
```

and so on, as according to {eq}`mmhoriz`.
```` 


The {term}`identity matrix` of size $n$, called $\mathbf{I}$ (or sometimes $\mathbf{I}_n$), is a diagonal $n\times n$ matrix with every diagonal entry equal to one. As can be seen from {eq}`mmhoriz` and {eq}`mmvert`, it satisfies $\mathbf{A}\mathbf{I}=\mathbf{A}$ for $\mathbf{A}\in\mathbb{C}^{m\times n}$ and $\mathbf{I}\mathbf{B}=\mathbf{B}$ for $\mathbf{B}\in\mathbb{C}^{n\times p}$. It is therefore the matrix analog of the number $1$.

````{prf:example}
Let
  
```{math}
  \mathbf{B} =
  \begin{bmatrix}
    2 & 1 & 7 & 4\\
    6 & 0 & -1 & 0\\
    -4 & -4 & 0 & 1
  \end{bmatrix}.
```

Suppose we want to create a zero in the (2,1) entry by adding $-3$ times the first row to the second row, leaving the other rows unchanged. We can express this operation as a product $\mathbf{A}\mathbf{B}$ as follows. From dimensional considerations alone, $\mathbf{A}$ will need to be $3\times 3$. According to {eq}`mvrow`, we get "$-3$ times row one plus row two" from left-multiplying $\mathbf{B}$ by the vector $\begin{bmatrix} -3 & 1 & 0 \end{bmatrix}$. Equation {eq}`mmvert` tells us that this must be the second row of $\mathbf{A}$. Since the first and third rows of $\mathbf{A}\mathbf{B}$ are the same as those of $\mathbf{B}$, similar logic tells us that the first and third rows of $\mathbf{A}$ are the same as the identity matrix:
  
```{math}
  \begin{bmatrix}
    1 & 0 & 0\\
    -3 & 1 & 0\\
    0 & 0 & 1
  \end{bmatrix} \mathbf{B}  =
  \begin{bmatrix}
    2 & 1 & 7 & 4\\
    0 & -3 & -22 & -12\\
    -4 & -4 & 0 & 1
  \end{bmatrix}.
```

This can be verified directly using {eq}`scalarmatrixmult`.
````

Note that a square matrix $\mathbf{A}$ can always be multiplied by itself to get a matrix of the same size. Hence we can define the integer powers $\mathbf{A}^2=(\mathbf{A})(\mathbf{A})$, $\mathbf{A}^3=(\mathbf{A}^2) \mathbf{A} = (\mathbf{A}) \mathbf{A}^2$ (by associativity), and so on. By definition, $\mathbf{A}^0=\mathbf{I}$.

According to the rules for multiplying matrices, there are two ways for vectors to be multiplied together. If $\mathbf{v}$ and $\mathbf{w}$ are in $\mathbb{C}^n$, their **inner product** is

```{math}
  \mathbf{v}^* \mathbf{w} = \sum_{k=1}^n \overline{v_k}\, w_k.
```

Trivially, one finds that $\mathbf{w}^* \mathbf{v} = \overline{(\mathbf{v}^*\mathbf{w})}$. Additionally, *any* two vectors $\mathbf{v}\in\mathbb{C}^m$ and $\mathbf{w}\in\mathbb{C}^n$ have an **outer product**, which is an $m\times n$ matrix:

```{math}
  \mathbf{v} \mathbf{w}^*
  =
  \begin{bmatrix}
    v_1\,\overline{w_1} & v_1\,\overline{w_2} & \cdots & v_1\,\overline{w_n}\\
    v_2\,\overline{w_1} & v_2\,\overline{w_2} & \cdots & v_2\,\overline{w_n}\\
    \vdots & \vdots & & \vdots\\
    v_m\,\overline{w_1} & v_m\,\overline{w_2} & \cdots & v_m\,\overline{w_n}
  \end{bmatrix}.
```

## Linear systems and inverses

Given a square, $n\times n$ matrix $\mathbf{A}$ and  $n$-vectors $\mathbf{x}$ and $\mathbf{b}$, the equation $\mathbf{A}\mathbf{x}=\mathbf{b}$ is equivalent to

```{math}
\begin{split}
  a_{11}x_1 + a_{12}x_2 + \cdots + a_{1n}x_n &= b_1 \\
  a_{21}x_1 + a_{22}x_2 + \cdots + a_{2n}x_n &= b_2 \\
  \vdots  \\
  a_{n1}x_1 + a_{n2}x_2 + \cdots + a_{nn}x_n &= b_n.
\end{split}
```

We say that $\mathbf{A}$ is {term}`invertible` or **nonsingular** if there exists another $n\times n$ matrix $\mathbf{A}^{-1}$, the {term}`inverse` of $\mathbf{A}$, such that $\mathbf{A}\mathbf{A}^{-1}=\mathbf{A}^{-1}\mathbf{A}=\mathbf{I}$, the identity matrix. Otherwise, $\mathbf{A}$ is **singular**. If a matrix is invertible, its inverse is unique. This and the following facts are usually proved in an elementary text on linear algebra.

(theorem-singularity)=

````{proof:theorem} (Linear algebra equivalence)
  The following statements are equivalent:

  1. $\mathbf{A}$ is nonsingular.
  2. $(\mathbf{A}^{-1})^{-1} = \mathbf{A}$.
  3. $\mathbf{A}\mathbf{x}=\boldsymbol{0}$ implies that $\mathbf{x}=\boldsymbol{0}$.
  4. $\mathbf{A}\mathbf{x}=\mathbf{b}$ has a unique solution, $\mathbf{x}=\mathbf{A}^{-1}\mathbf{b}$, for any $n$-vector $\mathbf{b}$.
````

## Exercises

1. ✍ In racquetball, the winner of a rally serves the next rally. Generally, the server has an advantage. Suppose that when Ashley and Barbara are playing racquetball, Ashley wins 60\% of the rallies she   serves and Barbara wins 70\% of the rallies she serves. If $\mathbf{x}\in\mathbb{R}^2$ is such that $x_1$ is the probability that Ashley serves first and $x_2=1-x_1$ is the probability that Barbara serves first, define a matrix $\mathbf{A}$ such that $\mathbf{A}\mathbf{x}$ is a vector of the probabilities that Ashley and Barbara each serve the second rally. What is the meaning of $A^{10}\mathbf{x}$?

2. ✍ Suppose we have lists of $n$ terms and $m$ documents. We can   define an $m\times n$ matrix $\mathbf{A}$ such that $A_{ij}=1$ if term $j$   appears in document $i$, and $A_{ij}=0$ otherwise. Now suppose that the term list is

    ``` julia
    "numerical", "analysis", "more", "cool", "accounting"
    ```

    and that $\mathbf{x} = \begin{bmatrix} 1 & 1 & 0 & 1 & 0  \end{bmatrix}^T$. Give an interpretation of the product $\mathbf{A}\mathbf{x}$.

    ````{only} solutions
    Entry $i$ of $A\mathbf{x}$ is the sum of $a_{ij}x_j$, so the result is the number of the terms 'numerical', 'analysis', and 'cool' that appear in document $i$. (Note: this is not the same as the  total number of occurrences of the words; each term is counted only once.)
    ````

3. ✍ Let
  
    ```{math}
    \mathbf{A} =
    \begin{bmatrix}
      0 & 1 & 0 & 0 \\
      0 & 0 & 0 & 1 \\
      0 & 0 & 0 & 0 \\
      0 & 0 & 1 & 0
    \end{bmatrix}.
    ```

    Show that $\mathbf{A}^n=0$ when $n\ge 4$.

4. ✍  Find two matrices $\mathbf{A}$ and $\mathbf{B}$, neither of which is the zero matrix, such that $\mathbf{A}\mathbf{B}=\boldsymbol{0}$.

   (prob-linalg-transposeidentity)=

5. ✍ Prove that when $\mathbf{A} \mathbf{B}$ is   defined, $\mathbf{B}^T\mathbf{A}^T$ is defined too, and use   equation {eq}`scalarmatrixmult` to show that   $(\mathbf{A}\mathbf{B})^T=\mathbf{B}^T\mathbf{A}^T$.

   (prob-linalg-inversetranspose)=

6. ✍ Show that if $\mathbf{A}$ is invertible, then  $(\mathbf{A}^T)^{-1}=(\mathbf{A}^{-1})^T$. (This matrix is often just written as $\mathbf{A}^{-T}$.)

7. ✍ Prove true, or give a counterexample: The product of upper triangular square matrices is upper triangular.
