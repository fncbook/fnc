# LU factorization

This is an alternate derivation of the LU factorization. It makes heavy use of the **outer product** of vectors.

## Outer products

Suppose $\mathbf{u}\in\real^m$ and $\mathbf{v}\in\real^n$. Then we define the outer product of these vectors as

```{math}
:label: outerproduct
\mathbf{u} \mathbf{v}^T =
\begin{bmatrix}
\mathbf{u}_1 v_1 & \mathbf{u}_1 v_2 & \cdots \mathbf{u}_1 v_n \\\mathbf{u}_2 v_1 & \mathbf{u}_2 v_2 & \cdots \mathbf{u}_2 v_n \\ \vdots & \vdots & & \vdots \\ \mathbf{u}_m v_1 & \mathbf{u}_m v_2 & \cdots \mathbf{u}_m v_n
\end{bmatrix}
```

There is an elegant way to express matrix multiplication entirely in terms of outer products. We illustrate here for a small example of $\mathbf{A}\in\rmn{2}{3}$, $\mathbf{B}\in\rmn{3}{3}$:

```{math}
\begin{split}
\begin{bmatrix} \mathbf{A}_{11} B_{11} + \mathbf{A}_{12}B_{21} + \mathbf{A}_{13}B_{31} & \mathbf{A}_{11} B_{12} + \mathbf{A}_{12}B_{22} + \mathbf{A}_{13}B_{32} & \mathbf{A}_{11} B_{13} + \mathbf{A}_{12}B_{23} + \mathbf{A}_{13}B_{33} \\ \mathbf{A}_{21} B_{11} + \mathbf{A}_{22}B_{21} + \mathbf{A}_{23}B_{31} & \mathbf{A}_{21} B_{12} + \mathbf{A}_{22}B_{22} + \mathbf{A}_{23}B_{32} & \mathbf{A}_{21} B_{13} + \mathbf{A}_{22}B_{23} + \mathbf{A}_{23}B_{33} \end{bmatrix}
& =
\begin{bmatrix} \mathbf{A}_{11} B_{11}  & \mathbf{A}_{11} B_{12}  & \mathbf{A}_{11} B_{13}  \\ \mathbf{A}_{21} B_{11} & \mathbf{A}_{21} B_{12} & \mathbf{A}_{21} B_{13}  \end{bmatrix}
+ \begin{bmatrix}  \mathbf{A}_{12}B_{21} &  \mathbf{A}_{12}B_{22}  &  \mathbf{A}_{12}B_{23}  \\  \mathbf{A}_{22}B_{21} &  \mathbf{A}_{22}B_{22}  &  \mathbf{A}_{22}B_{23}  \end{bmatrix}
+ \begin{bmatrix}  \mathbf{A}_{13}B_{31} & \mathbf{A}_{13}B_{32} & \mathbf{A}_{13}B_{33} \\ \mathbf{A}_{23}B_{31} & \mathbf{A}_{23}B_{32} & \mathbf{A}_{23}B_{33} \end{bmatrix} \\
& =
\begin{bmatrix} \mathbf{A}_{11} \\ \mathbf{A}_{21} \end{bmatrix}  \begin{bmatrix} B_{11} & B_{12} & B_{13} \end{bmatrix}
+ \begin{bmatrix} \mathbf{A}_{12} \\ \mathbf{A}_{22} \end{bmatrix}  \begin{bmatrix} B_{21} & B_{22} & B_{23} \end{bmatrix}
+ \begin{bmatrix} \mathbf{A}_{13} \\ \mathbf{A}_{23} \end{bmatrix}  \begin{bmatrix} B_{31} & B_{32} & B_{33} \end{bmatrix}
\end{split}
```

It is not hard to prove the generalization of this example to all matrix products. If we write the *columns* of $\mathbf{A}$ as $\mathbf{a}_1,\dots,\mathbf{a}_n$ and the *rows* of $\mathbf{B}$ as $\mathbf{b}_1^T,\dots,\mathbf{b}_n^T$, then

```{math}
:label: matrixouter
\mathbf{A}\mathbf{B} = \sum_{k=1}^n \mathbf{a}_k \mathbf{b}_k^T.
```

## Basic form

Our goal is to factor $\mathbf{A}=\mathbf{L}\mathbf{U}$, where $\mathbf{L}$ is unit lower triangular and $\mathbf{U}$ is upper triangular. Using the outer product form {eq}`matrixouter`, we have

```{math}
\mathbf{A} = \mathbf{L}\mathbf{U} = \sum_{k=1}^n \mathbf{\ell}_k \mathbf{u}_k^T.
```

The properties of $\mathbf{L}$ and $\mathbf{U}$ impose some useful structure:

```{math}
\begin{split}
\mathbf{A} = \begin{bmatrix} 1 \\ L_{21} \\ L{31} \\ \vdots \\ L_{n1} \end{bmatrix} \begin{bmatrix} \mathbf{u}_{11} & \mathbf{u}_{12} & \mathbf{u}_{13} & \vdots & \mathbf{u}_{1n} \end{bmatrix}
+ \begin{bmatrix} 0 \\ 1 \\ L_{32} \\ \vdots \\ L_{n2} \end{bmatrix} \begin{bmatrix} 0 & \mathbf{u}_{22} & \mathbf{u}_{23} & \vdots & \mathbf{u}_{2n} \end{bmatrix} \\
+ \begin{bmatrix} 0 \\ 0 \\ 1 \\ \vdots \\ L_{n3} \end{bmatrix} \begin{bmatrix} 0 & 0 & \mathbf{u}_{33} & \vdots & \mathbf{u}_{3n} \end{bmatrix}
+ \dots + \begin{bmatrix} 0 \\ 0 \\ 0 \\ \vdots \\ 1 \end{bmatrix} \begin{bmatrix} 0 & 0 & 0 & \vdots & \mathbf{u}_{nn} \end{bmatrix}.
\end{split}
```

All terms in the sum after the first are therefore zero in the first row and first column.

More compactly, we find

```{math}
\mathbf{e}_1^T \mathbf{A} = \sum_{k=1}^m (\mathbf{e}_1^T \mathbf{\ell}_k) \mathbf{u}_k^T = \mathbf{u}_1^T.
```

In words, the first row of $\mathbf{A}$ is equal to the first row of $\mathbf{U}$. Algorithmically this allows us to fill in that first row of $\mathbf{U}$. Similarly,

```{math}
\mathbf{A} \mathbf{e}_1 = \sum_{k=1}^m \mathbf{\ell}_k (\mathbf{u}_k^T \mathbf{e}_1) = \mathbf{u}_{11}\mathbf{\ell}_1.
```

Since $\mathbf{u}_{11}$ is already known, we now know how to find the first column of $\mathbf{L}$.

Now we make the definitions $\mathbf{A}_1=\mathbf{A}$ and

```{math}
\mathbf{A}_2 = \mathbf{A}_1 - \mathbf{\ell}_1 \mathbf{u}_1^T = \sum_{k=2}^m \mathbf{\ell}_k \mathbf{u}_k^T
```

We now play the same games with the second row and column of $\mathbf{A}_2$. Noting the property

```{math}
:label: Lcols
\mathbf{e}_j^T \mathbf{\ell}_{k} = L_{jk} = 0 \text{ if } j < k,
```

then for $j=2$ we find

```{math}
\mathbf{e}_2^T\mathbf{A}_2 = \sum_{k=2}^n (\mathbf{e}_2^T \mathbf{\ell}_k) \mathbf{u}_k^T = L_{22} \mathbf{u}_2^T = \mathbf{u}_2^T.
```

Therefore, the second row of $\mathbf{U}$ is the same as that of $\mathbf{A}_2$. We also apply

```{math}
:label: Urows
\mathbf{u}_{k}^T\mathbf{e}_j = \mathbf{u}_{kj} 0 \text{ if } j < k.
```

with $j=2$ to get

```{math}
\mathbf{A}_2 \mathbf{e}_2 = \sum_{k=2}^m \mathbf{\ell}_k (\mathbf{u}_k^T\mathbf{e}_2) = U_{22}\mathbf{\ell}_2.
```

We have just found $U_{22}$, so the second column of $\mathbf{A}_2$ is divided through to get $\mathbf{\ell}_2$.

Now we see how this is going to unfold. Starting at $j=1$, we calculate in order

```{math}
\mathbf{u}_j^T & = \mathbf{e}_j^T \mathbf{A}_j,\\
\mathbf{\ell}_j & = \frac{1}{U_{jj}} \mathbf{A}_j \mathbf{e}_j,\\
\mathbf{A}_{j+1} = \mathbf{A}_j - \mathbf{\ell}_j \mathbf{u}_j^T.
```

It's understood here that multiplication by a row or column of the identity matrix actually means to extract a row or column of $\mathbf{A}_j$.

INSERT FUNCTION HERE

## Pivoting

Above, we iteratively divided by the entries $U_{11},U_{22},\ldots$ as we found them. What if one of these were zero? There's a part of standard Gaussian elimination we have not yet used: swapping rows of the matrix.

We'll look at things a little differently. Rather than automatically choosing the first row and first column at the start, we will chose a different row, which we denote by $i_1$. To make things work, we need to replace the prior requirement in row 1 that

```{math}
L_{11} = 1, \quad L_{12}=L_{13}=\cdots=L_{1n}=0
```

with the equivalent in row $i_1$:

```{math}
L_{i_1,1} = 1, \quad L_{i_1,2}=L_{i_1,3}=\cdots=L_{i_1,n}=0.
```

Now our original first step of extracting the first row becomes

```{math}
\mathbf{e}_{i_1}^T \mathbf{A} = \sum_{k=1}^m (\mathbf{e}_{i_1}^T \mathbf{\ell}_k) \mathbf{u}_k^T = \mathbf{u}_1^T.
```

That is, the first row of $\mathbf{U}$ is equal to row $i_1$ of $\mathbf{A}_1$.

How do we decide how to select $i_1$? Now $U_{11}$ is the $(i_1,1)$ element of $\mathbf{A}_1$, which is known as the {term}`pivot element` in the first column. We need the pivot to be nonzero, since it will appear in a denominator. If no such element exists---i.e., if the entire first column is zero---then the matrix is singular, and the algorithm fails. For reasons we go into later, we select $i_1$ so that $|A_{i_1,1}|$ is maximized.  

Next, the first column extraction and derivation of $\mathbf{A}_2$ from $\mathbf{A}_1$ proceed as before. In $\mathbf{A}_2$ we find the row index $i_2$ that maximizes the absolute value over column 2, and row $i_2$ becomes row 2 of $\mathbf{U}$, etc. We update our algorithmic formulas to read

```{math}
\mathbf{u}_j^T & = \mathbf{e}_{i_j}^T \mathbf{A}_j,\\
\mathbf{\ell}_j & = \frac{1}{U_{jj}} \mathbf{A}_j \mathbf{e}_j,\\
\mathbf{A}_{j+1} = \mathbf{A}_j - \mathbf{\ell}_j \mathbf{u}_j^T.
```

That's easy enough to do, and we have still derived a factorization $\mathbf{A}=\mathbf{L}\mathbf{U}$. But the name $\mathbf{L}$ is no longer appropriate, because we destroyed the triangular structure when we used row $i_1$ in place of $1$. However, there is a closely related structure in its place.

To see it we need to look at an important detail left out up to this point: after $\mathbf{\ell}_1\mathbf{u}_1^T$ has been subtracted from $\mathbf{A}_1$, row $i_1$ and column 1 remain zero throughout the rest of the process. For example,

```{math}
\mathbf{A}_j \mathbf{e}_1 = \sum_{k=j}^m \mathbf{\ell}_k (\mathbf{u}_k^T \mathbf{e}_1) = \sum_{k=j}^m U_{k1} \mathbf{\ell}_k,
```

which is necessarily zero if $j>1$. A similar formula holds for rows. As a result, once $i_1$ is selected as the pivot row in column 1, it can never be selected again.

When we extend that reasoning throughout the process, we conclude that each row is selected exactly once as the pivot row. In other words, $i_1,\ldots,i_n$ is a **permutation** of the row indices $1,\dots,n$. And, if we shuffle the rows of $\mathbf{L}$ to put row $i_1$ first, then $i_2$, and so on, the result is once again a unit lower triangular matrix! We say that $\mathbf{L}$ is a **permuted** unit lower triangular matrix.

There's more. Suppose we rewind and apply the *inverse* permutation to the rows of $\mathbf{A}$: row 1 is moved to position $i_1$, row 2 is moved to row $i_2$, and so on. The algorithm will then play out exactly as before, but the pivot rows will be $i_1=1,i_2=2,\dots$. In that case $\mathbf{L}$ will again be truly lower triangular.

So we have derived a factorization of $\mathbf{A}$ into a permuted lower triangular and an upper triangular, or we equivalently have factored an inversely permuted $\mathbf{A}$ into two triangular factors. Algebraically these statements can be expressed as

```{math}
\begin{split}
\mathbf{A} &= \mathbf{P}^T \mathbf{L} \mathbf{U},\\
\mathbf{P} \mathbf{A} &= \mathbf{L} \mathbf{U},
```

where $\mathbf{P}$ is the {term}`permutation matrix` that results from shuffling the rows of the identity into positions $i_1,\dots,i_n$. A fact about permutation matrices that we shall not prove is that $\mathbf{P}^{-1} = \mathbf{P}^T$, which is why the two lines above are equivalent.

At the algorithmic level we need not bother with permutation matrices (and indeed it is less efficient to do so). 



which as before gives a way to extract $\mathbf{u}_1^T$. But now $ If we can't find an $i_1$ such that this is nonzero, then the entire first column of $A$ is zero, and this *would* imply that $A$ is singular. Otherwise, we have  exactly as before, and we know that we can compute $\mathbf{\ell}_1$. 

This is a lot less daunting than the formalism makes it sound. First, we use the maximum element in column 1 to select $i_1$ (more on this later).

```{code-cell}
i = zeros(Int,4);
i[1] = findmax(abs.(A1[:,1]))[2]
```

So we are targeting row 3 and column 1 to zero out.

```{code-cell}
A1 = A
U[1,:] = A1[i[1],:]
L[:,1] = A1[:,1]/U[1,1]
display(L)
display(U)
```

```{code-cell}
A2 = A1 - L[:,1]*U[1,:]'
```

Now we select a new row $i_2$ with a nonzero pivot in column 2.

```{code-cell}
i[2] = findmax(abs.(A2[:,2]))[2]
```

Now we want $\mathbf{e}_{i_2}^T\mathbf{A}_2=\mathbf{u}_2^T$. This happens if we require

$$L_{i_2,2} = 1, \quad L_{i_2,3}=\cdots=L_{i_2,m}=0.$$ 

(Note that $L_{i_2,1}$ was previously determined.)

```{code-cell}
U[2,:] = A2[i[2],:]
L[:,2] = A2[:,2]/U[2,2]
display(L)
display(U)
```

```{code-cell}
A3 = A2 - L[:,2]*U[2,:]'
```

By now the pattern is clear.

```{code-cell}
i[3] = findmax(abs.(A3[:,3]))[2]
U[3,:] = A3[i[3],:]
L[:,3] = A3[:,3]/U[3,3]
A4 = A3 - L[:,3]*U[3,:]'

i[4] = findmax(abs.(A4[:,4]))[2]
U[4,:] = A4[i[4],:]
L[:,4] = A4[:,4]/U[4,4];
```

Indeed, we did again get a factorization of $A$.

```{code-cell}
norm(A-L*U)
```

But what sort of factorization is it?

```{code-cell}
display(L)
display(U)
```

Just as before, $U$ is upper triangular. But $L$ is not triangular. However, think about the structural conditions imposed during the algorithm:

$$
L_{i_1,1} = 1, \quad L_{1_1,2}=\cdots=L_{i_1,m}=0,
$$ 

$$
L_{i_2,2} = 1, \quad L_{1_2,3}=\cdots=L_{i_2,m}=0,
$$ 

down to $L_{i_m,i_m}=1$ in the last step. What this means is that if we take the rows of $L$ in the order $i_1,i_2,i_3,i_4$, then the result is again unit lower triangular!

```{code-cell}
L[i,:]
```

We can express this result using a permutation matrix as $PL$. Conventionally though, this truly triangular matrix is the one we call $L$, and the one produced directly by the algorithm is $P^{-1}L=P^TL$. Since $A=P^TLU$, this implies that $PA=LU$. This is the **row-pivoted LU factorization** (or partially pivoted factorization, or $P^TLU$ factorization).

## Linear systems

The system $Ax=b$ is equivalent to $PAx=Pb$, or $L(Ux)=Pb$. We do a forward substitution using a permuted form of $b$, then a backward substitution using that result. (In practice we wouldn't move data around in memory, but just index the vector indirectly in the correct order.)

```{code-cell}
xact = ones(4);  b = A*xact;

Pb = similar(b);  Pb = b[i]; 
x = U\(L[i,:]\Pb)
```

The built-in `lu` returns a factorization object that will show you the (truly triangular) $L$ and $U$ factors, and the permutation vector we called `i` in the example above.

```{code-cell}
fact = lu(A);
typeof(fact)
```

```{code-cell}
fact.L
```

```{code-cell}
fact.U
```

```{code-cell}
fact.p
```
