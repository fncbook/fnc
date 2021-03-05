# Vector and matrix norms

The manipulations on matrices and vectors so far in this chapter have been algebraic, much like those in an introductory linear algebra course. In order to progress to the analysis of the algorithms we have introduced, we need a way to measure the "size" of vectors and matrices—size in the sense of magnitude or distance, not the number of rows and columns.

## Vector norms

```{index} norm; vector
```

For vectors we use a {term}`norm` $\| \cdot \|$, which is a function from $\real^n$ to $\real$ with the following properties:[^complexnorm]

[^complexnorm]: The same statements works for vectors with complex entries, with complex modulus in place of absolute values.
  
```{math}
:label: norm-properties
\begin{split}
\| \mathbf{x} \| &\ge 0 \text{ for all } \mathbf{x}\in\real^n \\
\| \mathbf{x} \| &=0 \text{ if and only if } \mathbf{x}=\boldsymbol{0} \\
\| \alpha \mathbf{x} \| &= |\alpha|\, \| \mathbf{x} \| \text{ for all } \mathbf{x},\,\alpha\\
\| \mathbf{x}+\mathbf{y} \| & \le \| \mathbf{x} \| + \| \mathbf{y} \| \text{ for all } \mathbf{x},\mathbf{y}
\end{split}
```

The last of these properties is known as the **triangle inequality**. Just as $\| \mathbf{x} \|=\| \mathbf{x}-\boldsymbol{0} \|$ is the distance from $\mathbf{x}$ to the origin, $\| \mathbf{x}-\mathbf{y} \|$ is the distance from $\mathbf{x}$ to $\mathbf{y}$. The three most important vector norms are

```{math}
\begin{split}
  \| \mathbf{x} \|_2 &= \left( \sum_{i=1}^n |x_i|^2 \right)^{\frac{1}{2}} = \sqrt{\mathbf{x}^T \mathbf{x}} \qquad \text{(2-norm)}   \\
  \| \mathbf{x} \|_\infty &= \max_{i=1,\dots,n} |x_i| \qquad \text{($\infty$-norm or max norm)}  \\
  \| \mathbf{x} \|_1 &= \sum_{i=1}^n |x_i| \qquad \text{(1-norm)}
\end{split}
```

```{proof:example} Julia demo
:class: demo
{doc}`demos/norms-vector`
```

The 2-norm corresponds to ordinary Euclidean distance. Most of the time, when just $\| \mathbf{x} \|$ is written, the 2-norm is implied. However, in this section we use it to mean a generic, unspecified vector norm.

```{index} unit vector
```

In any norm, we refer to a vector $\mathbf{x}$ satisfying $\| \mathbf{x} \|=1$ as a {term}`unit vector`. For any nonzero vector $\mathbf{v}$ we can find a unit vector through the normalization $\mathbf{x}=\mathbf{v}/\|\mathbf{v}\|$. Thus, we can interpret

```{math}
  :label: magdir
  \mathbf{v} = \| \mathbf{v} \| \: \frac{\mathbf{v}}{\| \mathbf{v} \|}.
```

as writing a nonzero vector $\mathbf{v}$ in magnitude--direction form. We say that a sequence of vectors $\mathbf{x}_n$ converges to $\mathbf{x}$ if

```{math}
  :label: vectorconverge
  \lim_{n\rightarrow\infty} \| \mathbf{x}_n - \mathbf{x} \| = 0.
```

By definition, a sequence is convergent in the infinity norm if and only if it converges componentwise. The same is true for a convergent sequence in *any* norm.

(theorem-normequivalence)=

```{proof:theorem} Norm equivalence
  In a finite-dimensional space, convergence in any norm implies convergence in all norms.
```

## Matrix norms

```{index} norm; matrix
```

```{index} norm; Frobenius
```

Although we view matrices as two-dimensional, we can also interpret them as vectors: simply stack the columns on top of one another.[^colmajor] Hence we can define matrix norms via vector norms. This is sometimes done with the vector 2-norm and leads to the matrix {term}`Frobenius norm`:

```{math}
:label: frobenius
\| \mathbf{A} \|_F = \left( \sum_{i,j} |A_{ij}|^2 \right)^{1/2}
```

This is the norm computed by the {term}`norm` function in Julia.

[^colmajor]: Column stacking is actually how matrices are stored in memory within Julia and is known as **column-major order**. MATLAB and Fortran also use column-major order, while C and Python use row-major order, in which the rows are stacked.

However, it often proves to be more useful to define matrix norms differently. Using a vector norm $\| \cdot \|_a$, we define for any $m\times n$ matrix $\mathbf{A}$,

```{math}
:label: matrixnorm
\| \mathbf{A} \|_{a} = \max_{\| \mathbf{x} \|_a=1} \| \mathbf{A}\mathbf{x} \|_a =
\max_{\mathbf{x}\neq \boldsymbol{0}} \frac{\| \mathbf{A}\mathbf{x} \|_a}{\| \mathbf{x} \|_a}.
```

(The last equality follows from linearity (as shown in [an exercise](problem-linearity).) This definition is called an {term}`induced matrix norm` and is computed by the {term}`opnorm` function in Julia. It is derived from the interpretation of a matrix as a linear operator between $\real^n$ and $\real^m$. Thus in the 2-norm, for instance,

```{math}
\| \mathbf{A} \|_2 = \max_{\| \mathbf{x} \|_2=1} \| \mathbf{A}\mathbf{x} \|_2.
```

For the rest of this section we will continue to omit subscripts when we want to refer to an unspecified norm; after this section, an unsubscripted norm is understood to be the 2-norm.

The definition of an induced matrix norm may seem overly complicated. However, there are some key properties that follow directly from the definition.

(theorem-norm-inequalities)=

````{proof:theorem} Norm inequalities
Let $\| \cdot \|$ designate a matrix norm and the vector norm that induced it. Then for all matrices and vectors of compatible sizes,

```{math}
:label: normineq1
\| \mathbf{A}\mathbf{x} \| \le \| \mathbf{A} \|\cdot \| \mathbf{x} \|.
```

For all matrices of compatible sizes,

```{math}
:label: normineq2
\| \mathbf{A}\mathbf{B} \| \le \| \mathbf{A} \|\cdot\| \mathbf{B} \|.
```

For a square matrix $\mathbf{A}$,

```{math}
:label: normineq3
\| \mathbf{A}^k \| \le \| \mathbf{A} \|^k \text{ for any integer $k\ge 0$.}
```
````

````{proof:proof}
The first result is trivial if $\mathbf{x}=\boldsymbol{0}$; otherwise,

```{math}
\frac{ \| \mathbf{A}\mathbf{x} \| }{\| \mathbf{x} \|} \le
\max_{\mathbf{x}\neq \boldsymbol{0}}  \frac{\| \mathbf{A}\mathbf{x} \|}{\| \mathbf{x} \|} = \| \mathbf{A} \|.
```

Inequality {eq}`normineq2` then follows because

```{math}
\| \mathbf{A}\mathbf{B}\mathbf{x} \| =\| \mathbf{A}(\mathbf{B}\mathbf{x}) \|\le \| \mathbf{A} \|\cdot\| \mathbf{B}\mathbf{x} \| \le
\| \mathbf{A} \|\cdot\| \mathbf{B} \|\cdot\| \mathbf{x} \|,
```

and then

```{math}
\| \mathbf{A}\mathbf{B} \| = \max_{\mathbf{x}\neq \boldsymbol{0}} \frac{\| \mathbf{A}\mathbf{B}\mathbf{x} \|}{\| \mathbf{x} \|} \le
\max_{\mathbf{x}\neq \boldsymbol{0}} \| \mathbf{A} \|\cdot\| \mathbf{B} \| = \| \mathbf{A} \|\cdot\| \mathbf{B} \|.
```

Finally,  {eq}`normineq3` results from repeated application of {eq}`normineq2`.
````

One can interpret the definition of an induced norm geometrically.  Each vector $\mathbf{x}$ on the unit "sphere" (as defined by the chosen vector norm) is mapped to its image $\mathbf{A}\mathbf{x}$, and the norm of $\mathbf{A}$ is the radius of the smallest "sphere" that encloses all such images. 

In addition, two of the vector norms we have encountered lead to equivalent formulas that are easy to compute from the matrix elements:

```{math}
:label: mxinfnorm
\| \mathbf{A} \|_\infty = \max_{1\le \,i \,\le n} \sum_{j=1}^n |A_{ij}|,
```

```{math}
:label: mxonenorm
\| \mathbf{A} \|_1 = \max_{1\le \,j\, \le n} \sum_{i=1}^n |A_{ij}|.
```

```{margin}
The infinity norm of a matrix is the maximum row sum, and the 1-norm is the maximum column sum, all after taking absolute values.
```

```{proof:example} Julia demo
:class: demo
{doc}`demos/norms-matrix`
```

In words, the infinity norm is the maximum *row* sum, and the 1-norm is the maximum *column* sum. (One way to keep this straight: the horizontal orientation of the $\infty$ symbol suggests the row sum for the $\infty$-norm, while the vertical orientation of the 1 suggests the column sum of the 1-norm.) In both cases you must take absolute values of the entries first! 

The geometric interpretation of the matrix 2-norm shown in {doc}`demos/norms-matrix`, as the radius of the smallest circle (or sphere or hypersphere in higher dimensions) containing the images of all unit vectors, is not a practical means of computing the norm. Nor is there a formula for it that makes it easy to compute directly from the matrix entries. The computation of the matrix 2-norm is discussed further in a later chapter.

## Exercises

1. ✍ Why is the vector 1-norm also called the "taxicab norm?"

2. ✍ **(a)** Draw the unit "circle" in the $\infty$-norm, i.e., the set of all vectors $\mathbf{x}\in\real^2$ such that $\| \mathbf{x} \|_\infty=1$.

    **(b)** Draw the unit "circle" in the 1-norm.
  
3. ✍ Show that for all vectors $\mathbf{x}\in\real^n$,

    **(a)** $\| \mathbf{x} \|_\infty \le \| \mathbf{x} \|_2$

    **(b)** $\| \mathbf{x} \|_2 \le \| \mathbf{x} \|_1$

    %%
    %% part (a)
    % Let $j$ be such that $|x_j|=\|x\|_\infty$. Then $\|x\|_2^2 = \sum x_i^2 \ge x_j^2 = \|x\|_\infty^2$.
    %% part (b)
    % $\|x\|_1^2 = \left( \sum_i |x_i| \right)^2 = \sum_{i,j} |x_i|\cdot |x_j| \ge \sum_i |x_i|^2 = \|x\|_2^2.$

4. ✍ Show that for any vectors $\mathbf{x}$, $\mathbf{y}$ in $\real^n$, $|\mathbf{x}^T\mathbf{y}| \le \| \mathbf{x} \|_1\| \mathbf{y} \|_\infty$.

    %%
    % 
    % $$|x^Ty| = \left| \displaystyle\sum_{i=1}^n x_i y_i \right| \le \displaystyle\sum_{i=1}^n |x_i| | y_i | $$
    % 
    % $$\le \left( \max_i |y_i| \right) \displaystyle\sum_{i=1}^n |x_i| =
    % \|y\|_\infty \|x\|_1$$

    (problem-linearity)=
5. ✍ Prove using the definition that for any induced matrix norm, matrix $\mathbf{A}$, and scalar $c$, $\| c\mathbf{A} \| = |c|\cdot \| \mathbf{A} \|$.

6. ✍ Let $\displaystyle \mathbf{A} =
    \begin{bmatrix}
      -1 & 1 \\ 2 & 2
    \end{bmatrix}$.

    **(a)** Find all vectors satisfying $\|\mathbf{x}\|_\infty=1$ and $\| \mathbf{A}\mathbf{x} \|_\infty=\| \mathbf{A} \|_\infty$.

    **(b)** Find a vector satisfying $\|\mathbf{x}\|_1=1$ and $\| \mathbf{A}\mathbf{x} \|_1=\| \mathbf{A} \|_1$.

    **(c)** Find a vector satisfying $\|\mathbf{x}\|_2=1$ such that $\| \mathbf{A}\mathbf{x} \|_2=\| \mathbf{A} \|_2$. (Hint: A unit 2-dimensional vector is a function only of its angle. Use the definition of $\|\mathbf{A}\|_2$ as the maximum of $\|\mathbf{A}\mathbf{x}\|_2$, which is a also a function of the angle.)

    %% (a)
    % The inf-norm of the matrix is 4. The relevant unit vectors are $[\pm 1;\pm1]$.
    %% (b)
    % The 1-norm of the matrix is 3. Some relevant unit vectors are $[\pm 1;0]$ and $[0;\pm 1]$, but there are others.
    %% (c)
    % If $x$ is $[\cos t;\sin t]$, then $Ax$ is $[\cos t-\sin t;2\cos t + 2\sin t]$, and $\|Ax\|_2^2$ is $1+4+6\cos t\sin t$. The maximum of this function over $t$ is 8, when e.g. $t=\pi/4$. So one of the answers is $x=[1/\sqrt{2};1/\sqrt{2}]$, and the other is its negative.

7. ✍ Prove the equivalence of the two formulas for a matrix norm in {eq}`matrixnorm`.

8. ✍ Explain why for any permutation matrix $\mathbf{P}$, $\| \mathbf{P} \|_2=1$.

    %%
    % The effect of $P$ is to rearrange the entries of any vector
    % $x$, so we conclude $\|Px\|_2=\|x\|_2$ for all $x$, which by
    % the definition of matrix norm implies $\|P\|_2=1$.    

9. ✍ Show that for any induced matrix norm and nonsingular matrix $\mathbf{A}$, $\| \mathbf{A}^{-1} \| \ge (\| \mathbf{A} \|)^{-1}$. (Hint: Apply [the norm inequalities theorem](theorem-norm-inequalities).)

10. ✍ **(a)** Show that for any $\mathbf{v}\in \real^n$,

    ```{math}
    \| \mathbf{v} \|_p \ge \max_{i=1,\ldots,n} |v_i|,
    ```

    where $p=1$, $2$, or $\infty$.

    **(b)** Show that for any $\mathbf{A}\in\real^{n \times n}$,

    ```{math}
    \| \mathbf{A} \|_p \ge \max_{i,j=1,\ldots,n} |A_{ij}|,
    ```

    where $p=1$, $2$, or $\infty$. (Hint: For $p=2$, rearrange {eq}`normineq1` for a well-chosen particular value of $\mathbf{x}$.)

    (problem-diagnorm)=
11. ✍ Show that if $\mathbf{D}$ is a diagonal matrix, then $\|\mathbf{D}\|_2 = \max_{i} |D_{ii}|$. You may assume the matrix is real and square, but that does not affect the result or the proof in any significant way. (Hint: Let $M=\max_{i} |D_{ii}|$. Proceed in two stages, showing that $\|\mathbf{D}\|_2\ge M$ and separately that $\|\mathbf{D}\|_2\le M$.)

    %%
    % Suppose that $x$ is any vector with $\|x\|_2=1$. The $i$ th entry of $Dx$
    % is $d_{ii}x_i$, so
    %
    % $$\left\|Dx\right\|^2 = \displaystyle\sum_i |d_{ii}|^2 \cdot
    % |x_i|^2  \le M^2  \displaystyle\sum_i |x_{i}|^2 = M^2,$$
    %
    % which proves that $\|D\|_2\le M$ by the definition of the 2-norm. 
    %
    % It remains to show that $\|D\|_2\ge M$. There must be some value of $k$
    % such that $|d_{kk}|=M$. Let
    % $$x=({\rm sign}\;d_{kk})\,e_k$.$
    % (Recall $e_k$ is column $k$ of the identity). Then
    % $Dx=({\rm sign}\;
    % d_{kk})^2 |d_{kk}|\, e_k = Me_k$ and $\|Dx\|_2=M$. Since $x$ is among the
    % unit vectors, maximizing over all unit vectors must be at least this
    % large; i.e., $\|D\|_2\ge M$.
    %
    % Hence $\|D\|_2=M$.

12. ✍ Suppose that $\mathbf{A}$ is ${n\times n}$ and that $\| \mathbf{A} \|<1$ in some induced matrix norm.

    **(a)** Show that $(\mathbf{I}-\mathbf{A})$ is nonsingular. (Hint: Show that $(\mathbf{I}-\mathbf{A})\mathbf{x}=\boldsymbol{0}$ for nonzero $\mathbf{x}$ implies that $\| \mathbf{A} \|\ge
    1$, using the definition of an induced matrix norm.)

    **(b)** Show that $\lim_{m\rightarrow \infty} \mathbf{A}^m = \boldsymbol{0}$. (For matrices as with vectors, we say $\mathbf{B}_m \rightarrow \mathbf{L}$ if $\| \mathbf{B}_m-\mathbf{L} \| \rightarrow 0$.)

    **(c)** Use (a) and (b) to show that we may obtain the geometric series

    ```{math}
    (\mathbf{I}-\mathbf{A})^{-1} = \sum_{k=0}^\infty \mathbf{A}^k.
    ```

    (Hint: Start with $\left(\sum_{k=0}^m \mathbf{A}^k\right)(\mathbf{I}-\mathbf{A})$ and take the limit.)
