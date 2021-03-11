# LU factorization

```{index} matrix factorization; LU
```

```{index} Gaussian elimination
```

Every first linear algebra course introduces {term}`Gaussian elimination` for a general square system $\mathbf{A}\mathbf{x}=\mathbf{b}$. In Gaussian elimination one uses row operations on an augmented matrix $[\mathbf{A}\; \mathbf{b}]$ to reduce it to an equivalent triangular system (usually upper triangular).

## Review example

Rather than writing out the process in full generality, we use an example to refresh your memory while getting arithmetic support from Julia.

```{prf:example} Julia demo
:class: demo
:label: demos-lu-gauss
{doc}`demos/lu-gauss`
```

## The algebra of Gaussian elimination

In [an earlier section](matrices) we observed that row and column operations can be expressed as linear algebra using columns from the identity matrix. This connection allows us to express Gaussian elimination using matrices. We will ignore the augmentation step, set aside $\mathbf{b}$ for now, and consider only the square system matrix $\mathbf{A}$.

As the first step in {prf:ref}`demos-lu-gauss`, we get the multiplier $A_{21}/A_{11}=-2$. The first row of $\mathbf{A}$ is extracted by $\mathbf{e}_1^T\mathbf{A}$. After $-2$ times this row is subtracted from row 2, with the other rows being left alone, we arrive at the matrix

```{math}
  \begin{bmatrix}
    \mathbf{e}_1^T \\[2mm] \mathbf{e}_2^T + 2 \mathbf{e}_1^T \\[2mm] \mathbf{e}_3^T \\[2mm] \mathbf{e}_4^T
  \end{bmatrix}\mathbf{A}.
```

This expression can be manipulated into

```{math}
\mathbf{A} =
  \left(\mathbf{I} +
  \begin{bmatrix}
    0 \mathbf{e}_1^T \\[2mm] 2 \mathbf{e}_1^T \\[2mm] 0 \mathbf{e}_1^T \\[2mm] 0 \mathbf{e}_1^T
  \end{bmatrix}
  \right) \mathbf{A}
  =\left(\mathbf{I} +2
  \begin{bmatrix}
    0 \\ 1 \\ 0 \\ 0
  \end{bmatrix}
  \mathbf{e}_1^T \right) \mathbf{A} = (\mathbf{I} + 2 \mathbf{e}_2\mathbf{e}_1^T) \mathbf{A}.
```

In general, adding $\alpha$ times row $j$ of $\mathbf{A}$ to row $i$ in place is done via the expression

```{math}
  :label: rowoperation
  (\mathbf{I} + \alpha \, \mathbf{e}_i \mathbf{e}_j^T ) \mathbf{A}.
```

Following many introductory texts on linear algebra, we refer to the matrix in parentheses above as an **elementary matrix**.

```{prf:example} Julia demo
:class: demo
:label: demos-lu-factors
{doc}`demos/lu-factors`
```

```{index} unit triangular matrix
```

The elementary matrix factors found in {prf:ref}`demos-lu-factors`, each in the form {eq}`rowoperation`, have some important properties. First, in addition to being triangular, each has all ones on the diagonal, so we call each a {term}`unit triangular matrix`. {prf:ref}`theorem-triangle-invert` implies that all unit triangular matrices are invertible, which is about to become important.

Let's review. The Gaussian elimination procedure in {prf:ref}`demos-lu-gauss` did six row operations in order to introduce six zeros into the lower triangle of $\mathbf{A}$. Each row operation can be expressed using multiplication by an elementary matrix $\mathbf{L}_{ij}$. At the end we get an upper triangular matrix, $\mathbf{U}$:

```{math}
:label: lufact0
\mathbf{U} = \mathbf{L}_{43}\mathbf{L}_{42}\mathbf{L}_{32}\mathbf{L}_{41}\mathbf{L}_{31}\mathbf{L}_{21} \mathbf{A}.
```

Now we multiply both sides on the left by $\mathbf{L}_{43}^{-1}$. On the right-hand side, it can be grouped together with $\mathbf{L}_{43}$ to form an identity matrix. Then we multiply both sides on the left by $\mathbf{L}_{42}^{-1}$, which knocks out the next term on the right side, etc. Eventually we get

```{math}
:label: lufact1
\mathbf{L}_{21}^{-1}\mathbf{L}_{31}^{-1}\mathbf{L}_{41}^{-1}\mathbf{L}_{32}^{-1}\mathbf{L}_{42}^{-1}\mathbf{L}_{43}^{-1} \mathbf{U} = \mathbf{A}.
```

We come next to an interesting property of these elementary matrices. If $i\ne j$, then for any scalar $\alpha$ we can calculate that

```{math}
\begin{split}
  \bigl(\mathbf{I} + \alpha\, \mathbf{e}_i \mathbf{e}_j^T\bigr) \bigl(\mathbf{I} - \alpha \,\mathbf{e}_i \mathbf{e}_j^T\bigr)
  & = \mathbf{I}  + \alpha \mathbf{e}_i \mathbf{e}_j^T - \alpha\, \mathbf{e}_i \mathbf{e}_j^T - \alpha^2\, \mathbf{e}_i \mathbf{e}_j^T \mathbf{e}_i \mathbf{e}_j^T \\
  &=  \mathbf{I}   - \alpha^2\, \mathbf{e}_i \bigl(\mathbf{e}_j^T \mathbf{e}_i\bigr)  \mathbf{e}_j^T = \mathbf{I},
\end{split}
```

since the inner product between any two different columns of $\mathbf{I}$ is zero. We have shown that

```{math}
   \bigl(\mathbf{I} + \alpha \, \mathbf{e}_i \mathbf{e}_j^T\bigr)^{-1} =  \mathbf{I} - \alpha\, \mathbf{e}_i \mathbf{e}_j^T.
```

All that is needed to invert $\mathbf{I} + \alpha\, \mathbf{e}_i \mathbf{e}_j^T$ is to flip the sign of $\alpha$, the lone element in the lower triangle.

We need one more remarkably convenient property of the elementary matrices. Looking, for example, at the first two matrices on the left in {eq}`lufact1`, we calculate that

```{math}
  (\mathbf{I} + \alpha \, \mathbf{e}_2\mathbf{e}_1^T) (\mathbf{I} + \beta\, \mathbf{e}_3\mathbf{e}_1^T) =  \mathbf{I} + \alpha\, \mathbf{e}_2\mathbf{e}_1^T + \beta \,\mathbf{e}_3\mathbf{e}_1^T + \alpha\beta \,\mathbf{e}_2\mathbf{e}_1^T\mathbf{e}_3\mathbf{e}_1^T.
```

We can use associativity in the last term to group together $\mathbf{e}_{1}^T\mathbf{e}_{3}$, which is another inner product that equals zero thanks to the structure of the identity matrix. In summary: the product of these elementary factors just combines the nonzero terms that each one has below the diagonal.

```{margin}
Gaussian elimination factors a square matrix into the product of a lower triangular matrix and an upper triangular matrix.
```

That reasoning carries across each of the new terms in the product on the left side of {eq}`lufact1`. The conclusion is that Gaussian elimination finds a unit lower triangular matrix $\mathbf{L}$ and an upper triangular matrix $\mathbf{U}$ such that

```{math}
:label: lufact
\mathbf{L} \mathbf{U} = \mathbf{A}.
```

Furthermore, the lower triangular entries of $\mathbf{L}$ are the row multipliers we found as in {prf:ref}`demos-lu-gauss`, and the entries of $\mathbf{U}$ are those found at the end of the elimination process. Equation {eq}`lufact` is called an {term}`LU factorization` of the matrix $\mathbf{A}$.

## An algorithm — for now

LU factorization reduces any linear system to two triangular ones. From this, solving $\mathbf{A}\mathbf{x}=\mathbf{b}$ follows immediately:

1. Factor $\mathbf{L}\mathbf{U}=\mathbf{A}$ using Gaussian elimination.
1. Solve $\mathbf{L}\mathbf{z}=\mathbf{b}$ for $\mathbf{z}$ using forward substitution.
1. Solve $\mathbf{U}\mathbf{x}=\mathbf{z}$ for $\mathbf{x}$ using backward substitution.

One of the important aspects of this algorithm is that the factorization step depends only on the matrix $\mathbf{A}$; the right-hand side $\mathbf{b}$ is not involved. Thus if one has to solve multiple systems with a single matrix $\mathbf{A}$, the factorization needs to be performed only once for all systems. As we show in [the next section](efficiency.md), the factorization is by far the most computationally expensive step, so this note is of more than academic interest.

Based on the examples and discussion above, a code for LU factorization is given in {numref}`Function {number}<function-lufact>`.

(function-lufact)=
````{proof:function} lufact
**LU factorization (not stable)**

```{code-block} julia
:lineno-start: 1
"""
lufact(A)

Compute the LU factorization of square matrix `A`, returning the
factors.
"""
function lufact(A)

n = size(A,1)
L = diagm(0=>ones(n))  # ones on main diagonal, zeros elsewhere
U = float(copy(A))

# Gaussian elimination
for j = 1:n-1
  for i = j+1:n
    L[i,j] = U[i,j] / U[j,j]   # row multiplier
    U[i,j:n] -= L[i,j]*U[j,j:n]
  end
end

return L,triu(U)
end
```
````

````{tip}
```{toggle}
Line 10 of {numref}`Function {number}<function-lufact>` points out two subtle Julia issues. First, arrays, including vectors and matrices, are really just references to blocks of memory. This data is much more efficient to pass around than the complete contents of the array. However, it means that a statement such as `U=A` just clones the array reference into the variable `U`. Any changes made to entries of `U` would then also be made to entries of `A`, because they refer to the same contents. In the context of `lufact` we don't want to change the original matrix, so we use {term}`copy` here to create an independent clone of the array contents and a new reference to them.

The second issue is that even when `A` has all integer entries, the LU factors may not. So if by copying `A` we create `U` as a matrix of integers, the function will fail with an {term}`InexactError` if we attempt to insert a noninteger result into `U` in line 16. To avoid this eventuality we convert `U` into a floating-point matrix with {term}`float`.
```
````

The multipliers are stored in the lower triangle of $\mathbf{L}$ as they are found. When operations are done to put zeros in column $j$, they are carried out only in lower rows to create an upper triangular matrix.  (Only columns $j$ through $n$ are accessed, since the other entries should already be zero.) At the end of the process the matrix $\mathbf{A}$ should be upper triangular, but since roundoff errors could create some small nonzeros the `triu` command is used to make them exactly zero.

```{prf:example} Julia demo
:class: demo
:label: demos-lu-function
{doc}`demos/lu-function`
```

Observe from {numref}`Function {number}<function-lufact>` that the factorization can fail if $A_{jj}=0$ when it is put in the denominator in line~14. This does *not* necessarily mean there is a zero in the diagonal of the original $\mathbf{A}$, because $\mathbf{A}$ is changed during the computation. Moreover, there are perfectly good nonsingular matrices for which this type of failure occurs, such as

```{math}
  \begin{bmatrix}
    0 & 1\\1 & 0
  \end{bmatrix}.
```

Fortunately, this defect can be repaired for all nonsingular matrices at minor cost, as we will see in [the section on pivoting](pivoting).

## Exercises

1. ✍ For each matrix, perform Gaussian elimination by hand to produce an LU factorization. Write out the $\mathbf{L}$ matrix using outer products of standard basis vectors.

    **(a)** $\displaystyle \begin{bmatrix}
    2 & 3 & 4 \\
    4 & 5 & 10 \\
    4 & 8 & 2
    \end{bmatrix}\qquad$
    **(b)** $\displaystyle \begin{bmatrix}
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

    **(b)** Find the LU factorization of $\mathbf{A}$.

    **(c)** Use the factors with triangular substitutions in order to solve $\mathbf{A}\mathbf{x}=\mathbf{b}$, and find $\mathbf{x}-\mathbf{z}$.
  
    ````{only} solutions
    ``` julia
    c = cos(pi/5);  s = sin(pi/5);
    A = [1 0 0; 0 1 0; 3 -1 1]*[c s 0;-s c 0;0 0 1]*[1 0 0;0 1 0;-3 1 1]
    z = [2, 2, 1];
    b = A*z
    [L,U] = FNC.lufact(A)
    x = FNC.backsub(U,FNC.forwardsub(L,b))
    x - z
    ```
    ````

    (problem-bigcorner)=
3. ⌨ In Julia, define
  
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

    **(a)** Using {numref}`Function {number}<function-lufact>` and triangular substitutions, solve the linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$, showing the result. About how many (to the nearest integer) accurate digits are in the result? (The answer is much less than the default 16 of double precision).

    **(b)** Repeat part (a) with $10^{20}$ as the corner element. (The result is even less accurate. We will study the causes of such low accuracy later in the chapter.)
  
    ````{only} solutions
    **(a)**
    ``` julia
    A = diag(ones(5,1)) + diag(ones(4,1),-1); A(1,5) = 1e12;
    x = (0:4)'/3;
    b = A*x;
    [L,U]=lufact(A);
    format long
    backsub(U,forwardsub(L,b))
    ```

    **(b)**
    ``` julia
    A = diag(ones(5,1)) + diag(ones(4,1),-1); A(1,5) = 1e20;
    x = (0:4)'/3;
    b = A*x;
    [L,U]=lufact(A);
    format long
    backsub(U,forwardsub(L,b))
    ```
    ````

4. not available
    % 4. \label{pro:dramadah} ⌨ Let $\mathbf{D}_n$ be the matrix created using MATLAB's {term}`gallery` function using where $n$ is a positive integer. It has interesting properties: the entries of $\mathbf{D}_n$ are all 0 or 1, and the entries of $\mathbf{D}_n^{-1}$ are all integers. Run an experiment that verifies that if $\mathbf{D}_n=\mathbf{L}\mathbf{U}$ is an LU factorization, then the entries of $\mathbf{L}$, $\mathbf{U}$, $\mathbf{L}^{-1}$, and $\mathbf{U}^{-1}$ are all integers for $n=2,3,\ldots,50$. You will have to do something more clever than visual inspection of the matrix entries to determine that they are integers; the {term}`round` and {term}`any` commands may be helpful.

    %for n = 2:50
    %    D = gallery('dramadah',n);
    %    [L,U] = lufact(D);
    %    Linv = inv(L);
    %    Uinv = inv(U);
    %    if any(any(L~=round(L))) || any(any(U~=round(U))) || ...
    %            any(any(Linv~=round(Linv))) || any(any(Uinv~=round(Uinv)))
    %        error
    %    end
    %end

5. ⌨ The {numref}`Function {number}<function-lufact>` function factors $\mathbf{A}=\mathbf{L}\mathbf{U}$ in such a way that $\mathbf{L}$ is a unit lower triangular matrix—that is, has all ones on the diagonal. It is also possible to define the factorization so that $\mathbf{U}$ is a unit upper triangular matrix instead. Write a function `lufact2` that uses {numref}`Function {number}<function-lufact>` *without modification* to produce this version of the factorization. (Hint: Begin with the standard LU factorization of $\mathbf{A}^T$.) Demonstrate on a nontrivial $4\times 4$ example.

    ````{only} solutions
    If $A=LU$, then $A^T=U^T L^T$, which is still of the form lower triangular $\times$ upper triangular. So if we apply our regular `lufact` to $A^T$, then the first matrix is $U^T$ and will leave $U$ as unit upper triangular.

    ``` julia
    A = round(8.7*magic(4));
    [L1,U1] = lufact(A');

    U = L1'

    L = U1'

    A - L*U
    ```

6. When computing the determinant of a matrix by hand, it's common to use cofactor expansion and apply the definition recursively. But this is terribly inefficient as a function of the matrix size.
  
    **(a)** ✍ Explain why, if $\mathbf{A}=\mathbf{L}\mathbf{U}$ is an LU factorization,

    ```{math}
      \det(\mathbf{A}) = U_{11}U_{22}\cdots U_{nn}=\prod_{i=1}^n U_{ii}.
    ```

    **(b)** ⌨ Using the result of part (a), write a function `determinant(A)` that computes the determinant using {numref}`Function {number}<function-lufact>`. Test your function on at least two nontriangular $5\times 5$ matrices, comparing your result to the result of `det` from the standard `LinearAlgebra` package.

    <!-- Use your function and the built-in "det" on the matrices "magic(n)" for $n=3,4,\ldots,7$, and make a table showing $n$, the value from your function, and the relative error when compared to "det". -->
  
    <!-- **(c)** ⌨ Show that "determinant" fails for "magic(8)" but is fine for "magic(9)". Speculate on what property of these two matrices makes the results so different. -->
  
7. not available

    <!-- ⌨ \label{pro:LUoneloop} Consider the portion of {numref}`Function {number}<function-lufact>` in the innermost
      loop. Because the different iterations in $i$ are all independent,
      it is possible to rewrite this group of operations without a
      loop. In fact, the necessary changes are to delete the keyword
      "for" in the inner loop, and delete the following
      "end" line. (You should also put a semicolon at the end of
      "i = j+1:n" to suppress extra output.)
      
      1. 
      \item Make the changes as directed and verify that the function
        works properly on the matrix from {prf:ref}`demos-lufactfun`.
      \item Write out symbolically (i.e., using ordinary elementwise
        vector and matrix notation) what the new version of the function
        does in the case $n=5$ for the iteration with $j=3$.
       -->
