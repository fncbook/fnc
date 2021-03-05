# Row pivoting

```{index} matrix factorization; LU
```

```{proof:example} Julia demo
:class: demo
{doc}`demos/pivot-fail`
```

As mentioned in [the section on LU factorization](lu), the $\mathbf{A}=\mathbf{L}\mathbf{U}$ factorization
does not work for every nonsingular $\mathbf{A}$. A simple rearrangement of the system found in {doc}`demos/lu-gauss` demonstrates this fact.

The breakdown in {doc}`demos/pivot-fail` is easily explained. After elimination is performed successfully in the first column, we have the matrix

```{math}
  \begin{bmatrix}
    2 & 0 & 4 & 3 \\
    0 & 0 & 6 & -10 \\
    0 & 15 & 0 & -6 \\
    0 & 5 & 1 & -4 \\
  \end{bmatrix}.
```

At this point there is no way to proceed, because trying to find the next row multiplier means division by the zero now sitting in the (2,2) position. In this particular case, we know that we can remove the problem by swapping two of the rows. Fortunately, that capability is all that is needed for *any* invertible matrix.

```{index} pivoting
```

A row swap is necessary when the algorithm would otherwise require a division by zero. The only context in which this can occur is if the $(j,j)$ element is zero just before the eliminations for column $j$. This number is called the **pivot element**. If it is zero, we need to swap row $j$ with another row in order to create a nonzero pivot. But we don't want to swap with a row *above* row $j$, because that could destroy the upper triangular structure we have partially built. So we search only in rows $j+1$ through $n$ of column $j$ to find a nonzero element to serve as pivot. This process is called {term}`row pivoting` or *partial pivoting*.

What if no nonzero pivot can be found? The following theorem ties up this loose end.

(theorem-pivot)=

```{proof:theorem} Row pivoting
  If a pivot element and all the elements below it are zero, then the original matrix $\mathbf{A}$ is singular. In other words, if $\mathbf{A}$ is nonsingular, then Gaussian elimination with row pivoting will run to completion.
```

```{proof:proof}
The proof of the first statement is considered in [an exercise](problem-zerocolumn). The second statement follows from it, because the only way that the Gaussian elimination algorithm can have an undefined step is when all of the choices for a  pivot element are zero.
```

It's important to keep in mind that [the pivoting theorem](theorem-pivot) does not refer to elements of the *original* matrix $\mathbf{A}$. Only the numbers that appear during the elimination process are relevant.

## The algebra of row pivoting

Recall that the algebraic expression of LU factorization culminated in equations {eq}`lufact0` and {eq}`lufact1`. Pivoting introduces a new type of elementary matrix called a {term}`permutation matrix`, which is an identity matrix with its rows (or depending on your point of view, its columns) reordered. The structure of the elimination process implies

```{math}
:label: pivot0
\mathbf{U} = \mathbf{L}_{43}\mathbf{P}_3(\mathbf{L}_{42}\mathbf{L}_{32})\mathbf{P}_2 (\mathbf{L}_{41}\mathbf{L}_{31}\mathbf{L}_{21})\mathbf{P}_1 \mathbf{A},
```

for permutation matrices $\mathbf{P}_j$ each featuring a single row swap. With a little more algebra that we won't give here, this equation can be rearranged to

```{math}
:label: PALU
\mathbf{P}\mathbf{A} = \mathbf{L}\mathbf{U},
```

where $\mathbf{P}$ is the permutation matrix resulting from all of the row swaps in the factorization, and $\mathbf{L}$ is unit lower triangular. (The construction of $\mathbf{L}$ also has to account for the row swaps; we do not give the details.) We emphasize that $\mathbf{P}$ cannot be computed all at once from inspection of the original elements of $\mathbf{A}$. Rather, the elimination process has to be run in full in order to deduce the row swaps.

```{index} matrix factorization; PLU
```

```{proof:example} Julia demo
:class: demo
{doc}`demos/pivot-fixed`
```

Permutation matrices have the property $\mathbf{P}^T = \mathbf{P}^{-1}$ (see [this exercise](problem-perminverse)). As a result, we can rewrite equation {eq}`PALU` as $\mathbf{A}=\mathbf{P}^T\mathbf{L}\mathbf{U}$, and we have a modified factorization of $\mathbf{A}$ (which we call a {term}`PLU factorization`). Either way, we could rewrite the linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$ as $\mathbf{L}\mathbf{U}\mathbf{x}=\mathbf{P}\mathbf{b}$, which implies that we can use backward and forward substitution following a permutation of the elements of $\mathbf{b}$. However, Julia offers an alternative shortcut for this situation.

Row pivoting does not add any flops to the factorization process. It does introduce $O(n^2)$ numerical comparisons, which should be inconsequential compared to the $O(n^3)$ flop requirement.

## Stability and pivoting

```{index} stability
```

If you're very attentive, you may have noticed something curious in the demo {doc}`demos/pivot-fixed`. The matrix we introduced had a single swap between rows 2 and 4, and therefore the elimination process could have been fixed by making just that swap again. But this is not what is shown by the $\mathbf{P}$ computed by Julia. The reason is that there is an additional rule about row pivoting that is essential to the success of the algorithm in the presence of rounding errors.

A simple example illustrates the key issue. Let

```{math}
  \mathbf{A} =
  \begin{bmatrix}
    -\epsilon & 1 \\ 1 & -1
  \end{bmatrix}, \qquad
  \mathbf{b} =
  \begin{bmatrix}
    1-\epsilon \\ 0
  \end{bmatrix},
```

where $\epsilon$ is a small positive number.  Elimination can be performed  without any pivoting:

```{math}
  \begin{bmatrix}
    -\epsilon & 1 & 1-\epsilon\\
    0 & -1+\epsilon^{-1} & \epsilon^{-1} - 1
  \end{bmatrix} \quad \Rightarrow \quad
  \begin{split}
    x_2 &= 1\\
    x_1 &= \frac{(1-\epsilon) - 1}{-\epsilon}.
  \end{split}
```

```{index} subtractive cancellation
```

In exact arithmetic, this produces the correct solution $x_1=x_2=1$.  But look at how $x_1$ is computed. It involves the subtraction of nearby numbers, which is sure to lead to subtractive cancellation in floating-point arithmetic. As a result, the solution will lose an arbitrary amount of accuracy, depending on the value of $\epsilon$. In the language of [the section on stability](../intro/stability), the computation has a poorly conditioned step, introducing instability into the complete algorithm.

Now suppose we swapped the rows of the matrix before elimination, even though it isn't actually required.

```{math}
  \begin{bmatrix}
    1 & -1 & 0\\
    0 & 1-\epsilon & 1-\epsilon
  \end{bmatrix} \quad \Rightarrow \quad
  \begin{aligned}
    x_2 &= 1\\
    x_1 &= \frac{0 - (-1)}{1}.
  \end{aligned}
```

Each of the arithmetic steps is well-conditioned now, and the solution is computed stably.

In general, a small (in absolute value) pivot element means weak dependence on the variable that is about to be eliminated. Even though only a zero pivot makes elimination technically impossible, in floating point using a pivot close to zero can cause instability. As a result, the pivoting rule chosen in practice is the following:

```{proof:myrule} Pivoting
When performing elimination in column $j$, swap row $j$ with the row below it whose entry in column $j$ is the largest (in absolute value).
```

One consequence of this rule is that all the below-diagonal elements in the unit lower triangular matrix $\mathbf{L}$ are bounded above by 1 in absolute value.

## Exercises

1. ✍ Let $\mathbf{A}$ be a square matrix and $\mathbf{b}$ be a column vector of compatible length. Here is correct Julia code to solve $\mathbf{A}\mathbf{x}=\mathbf{b}$:

    ``` julia
    using LinearAlgebra
    L,U = lu(A)
    x = U \ (L\b)
    ```

    Suppose instead you replace the last line above with

    ``` julia
    x = U \ L \ b
    ```

    Mathematically, what vector is found?

    %%
    % The difference between |U\(L\b)| and |U\L\b| has to do with the ordering
    % of evaluation. The first case does one forward substitution, then
    % one backward substitution, as desired. The second case does $n$ forward
    % substitutions using the columns of $L$ as right-hand sides. This matrix
    % does not even have any particular structure, for example,
    %[L,U] = lu(magic(5));
    %U \ L

    %%
    % So then the next step is to do a linear system solution, requiring
    % another _LU_ factorization, plus triangular solutions! Mathematically,
    % the backslash is always like multiplying on the left by the inverse, so
    % the operation performed is
    %
    % $$(U^{-1}L)^{-1}b=L^{-1}Ub,$$
    %
    % which is not the same as the correct solution, $U^{-1}L^{-1}b$.

2. ⌨ Rework {doc}`demos/lu-gauss` with row pivoting, using Julia to do the elimination step by step and visual inspection of the values as you go to determine the correct pivot element for each column. Start with $\mathbf{P}=\mathbf{I}$ and perform all the same row swaps in $\mathbf{P}$ as you do in the process of creating $\mathbf{U}$. (There is no need to compute $\mathbf{L}$.)  

    (problem-lumpstring)=
3. Suppose a string is stretched with tension $\tau$ horizontally between two anchors at $x=0$ and $x=1$. At each of the $n-1$ equally spaced positions $x_k=k/n$, $k=1,\ldots,n-1$, we attach a little mass $m_i$ and allow the string to come to equilibrium. This causes vertical displacement of the string. Let $q_k$ be the amount of displacement at $x_k$. If the displacements are not too large, then an approximate force balance equation is
  
    ```{math} n \tau (q_k - q_{k-1}) + n\tau (q_k - q_{k+1}) =
    m_k g, \qquad k=1,\ldots,n-1,
    ```

    where $g=-9.8$ m/s$^2$ is the acceleration due to gravity, and we define $q_0=0$ and $q_n=0$ due to the anchors.

    **(a)** ✍ Show that the force balance equations can be written as a linear system $\mathbf{A}\mathbf{q}=\mathbf{f}$, where $\mathbf{q}$ is a vector of displacements and $\mathbf{A}$ is a tridiagonal matrix, i.e., $A_{ij}=0$ if $|i-j|>1$) of size $(n-1)\times(n-1)$.

    **(b)** ⌨  Let $\tau=10$ N, and $m_k=(1/10n)$ kg for every $k$. Find the displacements when $n=4$ and $n=40$, and superimpose plots of $\mathbf{q}$ over $0\le x \le 1$ for the two cases. (Be sure to include the zero values at $x=0$ and $x=1$ in your plots of the string.)

    **(c)** ⌨  Repeat (b) for the case $m_k = (k/5n^2)$ kg.

    %% part (a)
    % The unknowns in the system are $q_1,\ldots,q_{n-1}$. Entry $a_{ij}$ is
    % $2n\tau$ if $i=j$ and $-n\tau$ if $|i-j|=1$. Entry $i$ of the right-side 
    % vector is $gm_i$.

    %% part(b)
    <!-- tau = 10;
    for n = [4,40]
    m = 1/(10*n)*ones(n-1,1);
    f = -9.8*m;
    A = toeplitz([2*n*tau,-n*tau,zeros(1,n-3)]);
    q = A\f;
    x = (0:n)'/n;
    subplot(1,2,1+floor(n/40)), plot(x,[0;q;0],'-o')
    end -->

    %% part (c)
    <!-- tau = 10;
    for n = [4,40]
    m = (1:n-1)' ./ (5*n.^2);
    f = -9.8*m;
    A = toeplitz([2*n*tau,-n*tau,zeros(1,n-3)]);
    q = A\f;
    x = (0:n)'/n;
    subplot(1,2,1+floor(n/40)), plot(x,[0;q;0],'-o')
    end -->

4. not available

    %⌨ Repeat the experiment of {ref}`prob-lufact-dramadah` using $\mathbf{P}\mathbf{A}=\mathbf{L}\mathbf{U}$ in place of $\mathbf{A}=\mathbf{L}\mathbf{U}$. Does the hypothesis on the entries of $\mathbf{L}$, $\mathbf{U}$, and their inverses remain true? 

    (problem-zerocolumn)=
5. ✍ Suppose that $\mathbf{A}$ is an $n\times n$ matrix such that its first $k$ columns are upper triangular and $A_{kk}=0$. Show that $\mathbf{A}$ is singular. (This can be used to complete the proof of [the pivoting theorem](theorem-pivot), because such a matrix occurs when there is no possible choice of a nonzero pivot.)

6. ✍ How many unique $6\times 6$ permutation matrices are there?

7. ✍ Suppose that `A` is a $4\times 6$ matrix in Julia and you define

    ``` julia
    B = A[end:-1:1,end:-1:1]
    ```

    Show that $\mathbf{B} = \mathbf{P} \mathbf{A} \mathbf{Q}$ for certain matrices $\mathbf{P}$ and $\mathbf{Q}$.

    (problem-perminverse)=
8. ✍ Suppose an $n\times n$ permutation matrix $\mathbf{P}$ has the effect of moving rows $1,2,\ldots,n$ to new positions $i_1,i_2,\ldots,i_n$. Then $\mathbf{P}$ can be expressed as
  
    ```{math}
    \mathbf{P} = \mathbf{e}_{i_1}\mathbf{e}_1^T + \mathbf{e}_{i_2}\mathbf{e}_2^T + \cdots + \mathbf{e}_{i_n}\mathbf{e}_n^T.
    ```

    **(a)** For the case $n=4$ and $i_1=3$, $i_2=2$, $i_3=4$, $i_4=1$, write out separately, as matrices, all four of the terms in the sum. Then add them together to find $\mathbf{P}$.

    **(b)** Use the formula in the general case to show that $\mathbf{P}^{-1}=\mathbf{P}^T$.
  
    %%
    %% (a)
    <!-- I = eye(4);
    P1 = I(:,3)*I(:,1)'
    P2 = I(:,2)*I(:,2)'
    P3 = I(:,4)*I(:,3)'
    P4 = I(:,1)*I(:,4)'
    P = P1 + P2 + P3 + P4 -->
    %% (b)
    % Form $P^T P$. Writing out all of the terms, associate the two inner
    % factors to get an inner product. This is zero in every case except
    % when $i_k=k$, when it is one. So you add up $e_ke_k^T$ over all $k$,
    % which is the identity.  
