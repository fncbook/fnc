# Computing QR factorizations

```{index} orthogonal; matrix
```

To compute an LU factorization, we follow elimination rules to introduce zeros into the lower triangle of the matrix, leaving only the $\mathbf{U}$ factor. The row operations are themselves triangular and can be combined into the $\mathbf{L}$ factor. For the QR factorization, the game is again to introduce zeros into the lower triangle, but the rules have changed; now the row operations must be done orthogonally. Thanks to [the orthogonal matrix theorem](theorem-orthogmatrix), the product of orthogonal operations will also be orthogonal, providing us with the final $\mathbf{Q}$ of the QR.

## Householder reflections

```{index} Householder reflector
```

```{index} orthogonal; matrix
```

A **Householder reflector** is a particular type of orthogonal matrix $\mathbf{P}$. The reflection is customized for a particular given vector $\mathbf{z}$ so that $\mathbf{P}\mathbf{z}$ is nonzero only in the first element.  Since orthogonal matrices preserve the 2-norm, we must have

```{math}
:label: hhgoal
\mathbf{P}\mathbf{z} =
\begin{bmatrix}
\pm \| \mathbf{z} \|\\0 \\ \vdots \\ 0
\end{bmatrix} = \pm \| \mathbf{z} \| \mathbf{e}_1.
```

(Recall that $\mathbf{e}_k$ is the $k$th column of the identity matrix.) We choose the positive sign for our discussion, but see [`qrfact`](function-qrfact) and [an exercise](problem-hhsign) for important computational details.

Given $\mathbf{z}$, let

```{math}
  :label: hhvector
  \mathbf{v} = \| \mathbf{z} \| \mathbf{e}_1-\mathbf{z}.
```

Then the reflector $\mathbf{P}$ is defined by

```{math}
  :label: hhreflect
  \mathbf{P} = \mathbf{I} - 2\frac{\mathbf{v} \mathbf{v}^T}{\mathbf{v}^T\mathbf{v}}
```

if $\mathbf{v}\neq\boldsymbol{0}$, or $\mathbf{P}=\mathbf{I}$ if $\mathbf{v}=\boldsymbol{0}$. Note that $\mathbf{v}^T\mathbf{v}$ is a scalar and can appear in a denominator, while the outer product $\mathbf{v} \mathbf{v}^T$ is $n\times n$. It is straightforward to show that $\mathbf{P}$ has the following key properties.

```{index} symmetric matrix
```

(theorem-hhreflect)=

````{proof:theorem} (Householder reflector)
Let $\mathbf{v}=\| \mathbf{z} \|\mathbf{e}_1-\mathbf{z}$ and let $\mathbf{P}$ be given by {eq}`hhreflect`. Then $\mathbf{P}$ is symmetric and orthogonal, and $\mathbf{P}\mathbf{z}=\| \mathbf{z} \|\mathbf{e}_1$.
````

````{proof:proof}
The case $\mathbf{v}=\boldsymbol{0}$ is obvious. For $\mathbf{v}\neq\boldsymbol{0}$, the proofs of symmetry and orthogonality are left to the exercises. As for the last fact, we simply compute
  
```{math}
:label: hhapply
\mathbf{P}\mathbf{z} = \mathbf{z} - 2 \frac{\mathbf{v} \mathbf{v}^T \mathbf{z}}{\mathbf{v}^T\mathbf{v}}
= \mathbf{z} - 2 \frac{\mathbf{v}^T \mathbf{z}}{\mathbf{v}^T\mathbf{v}} \mathbf{v},
```

and, since $\mathbf{e}_1^T\mathbf{z}=z_1$,
  
```{math}
\begin{split}
    \mathbf{v}^T\mathbf{v} &= \| \mathbf{z} \|^2 - 2 \| \mathbf{z} \| z_1 + \mathbf{z}^T\mathbf{z}
    = 2\| \mathbf{z} \|(\| \mathbf{z} \|-z_1),\\
    \mathbf{v}^T\mathbf{z} &= \| \mathbf{z} \|z_1 - \mathbf{z}^T\mathbf{z} = -\| \mathbf{z} \|\bigl(\| \mathbf{z} \|-z_1\bigr),
\end{split}
```

leading finally to
  
```{math}
\mathbf{P}\mathbf{z} = \mathbf{z} - 2\cdot
\frac{-\| \mathbf{z} \| \bigl(\| \mathbf{z} \|-z_1\bigr)}{2\| \mathbf{z} \| \bigl(\| \mathbf{z} \|-z_1\bigr)} \mathbf{v}
= \mathbf{z} + \mathbf{v} = \| \mathbf{z} \|\mathbf{e}_1.
```
````

The reason $\mathbf{P}$ is called a reflector is sketched in {numref}`fig-hhreflect`.

```{figure} figures/hhreflect.svg
:name: fig-hhreflect
Householder reflector.
```

```{index} orthogonal
```

The vector $\mathbf{v}$ defines an $n-1$ dimensional subspace $S$ perpendicular to it. Elementary vector analysis shows that $\mathbf{v}(\mathbf{v}^T\mathbf{z})/(\mathbf{v}^T\mathbf{v})$ is the vector projection of $\mathbf{z}$ along $\mathbf{v}$. By subtracting this projection from $\mathbf{z}$, we end up in $S$, but by subtracting twice the projection we get a reflection through $S$. This reflection occurs when $\mathbf{P}$ is applied to any vector, but when it is applied to $\mathbf{z}$, the result ends up on the $x_1$-axis.

## Factorization algorithm

```{index} matrix factorization; QR
```

```{margin}
The QR factorization is computed by using successive Householder reflections to introduce zeros in one column at a time.
```

````{proof:example} Julia demo
:class: demo
{doc}`demos/house-qr`
````

The QR factorization is computed by using successive Householder reflections to introduce zeros in one column at a time. We first show the process for a small numerical example in {doc}`demos/house-qr`.

You may be wondering what happened to the $\mathbf{Q}$ in {doc}`demos/house-qr`. Each Householder reflector is orthogonal but not full-size. We have to pad it out to represent algebraically the fact that a block of the first rows are left alone. Given a reflector $\mathbf{P}_k$ that is of square size $m-k+1$, we define

```{math}
\mathbf{Q}_k =
\begin{bmatrix}
\mathbf{I}_{k-1} & \boldsymbol{0} \\ \boldsymbol{0} & \mathbf{P}_k
\end{bmatrix}.
```

It is trivial to show that $\mathbf{Q}_k$ is also orthogonal. Then

```{math}
  \mathbf{Q}_n \mathbf{Q}_{n-1}\cdots \mathbf{Q}_1 \mathbf{A} = \mathbf{R}.
```

But $\mathbf{Q}_n \mathbf{Q}_{n-1}\cdots \mathbf{Q}_1$ is orthogonal too, and we multiply on the left by its transpose to get $\mathbf{A}=\mathbf{Q}\mathbf{R}$, where $\mathbf{Q} =  (\mathbf{Q}_n \mathbf{Q}_{n-1}\cdots \mathbf{Q}_1)^T$. We don't even need to form these matrices explicitly. Writing

```{math}
  \mathbf{Q}^T = \mathbf{Q}_n \mathbf{Q}_{n-1}\cdots \mathbf{Q}_1 = \mathbf{Q}_n \Bigl( \mathbf{Q}_{n-1}\bigl(\cdots (\mathbf{Q}_1\mathbf{I})\cdots\bigr)\Bigr),
```

we can build $\mathbf{Q}^T$ iteratively by starting with the identity and doing the same row operations as on $\mathbf{A}$.  Creating $\mathbf{Q}^T$ with row operations on $\mathbf{I}$ uses much less memory than building the $\mathbf{Q}_k$ matrices explicitly.

The algorithm we have described is encapsulated in [`qrfact`](function-qrfact). There is one more refinement in it, though. As indicated by {eq}`hhapply`, the application of a reflector $\mathbf{P}$ to a vector does not require the formation of the matrix $\mathbf{P}$ itself. Its effect can be computed directly from the vector $\mathbf{v}$, as is shown in [`qrfact`](function-qrfact).

(function-qrfact)=

```{proof:function} qrfact

**QR factorization by Householder reflections.**

```{code-block} julia
:lineno-start: 1
"""
qrfact(A)

QR factorization by Householder reflections. Returns Q and R.
"""
function qrfact(A)

    m,n = size(A)
    Qt = Matrix(Diagonal(ones(m)))
    R = float(copy(A))
    for k in 1:n
      z = R[k:m,k]
      v = [ -sign(z[1])*norm(z) - z[1]; -z[2:end] ]
      nrmv = norm(v)
      if nrmv < eps() continue; end  # skip this iteration
      v = v / nrmv;                  # simplifies other formulas
      # Apply the reflection to each relevant column of A and Q
      for j in 1:n
        R[k:m,j] -= v*( 2*(v'*R[k:m,j]) )
      end
      for j in 1:m
        Qt[k:m,j] -= v*( 2*(v'*Qt[k:m,j]) )
      end
    end
    
    return Qt',triu(R)
end
```

The Julia command `qr` works similarly to, but more efficiently than, [`qrfact`](function-qrfact). It finds the factorization in $\sim(2m n^2-n^3/3)$ flops asymptotically.

## Exercises

1. ⌨ Find a Householder reflector $\mathbf{P}$ such that
  
    ```{math}
    \mathbf{P}
    \begin{bmatrix}
      -6 \\ 2 \\ 9
    \end{bmatrix} =
    \begin{bmatrix}
      11\\0\\0
    \end{bmatrix}.
    ```

2. ✍ Prove the unfinished items in the [reflector theorem](theorem-hhreflect), namely that a Householder reflector $\mathbf{P}$ is symmetric and orthogonal.

    % Prove that the Householder reflector $P$ is symmetric and orthogonal.  First,
    % we need to show that $P^T=P$, using
    % \[ P = I - \frac{2}{v^Tv} vv^T,\]
    % we have
    % \begin{eqnarray*}
    % P^T & = &  \left( I - \frac{2}{v^Tv} vv^T \right)^T \\
    % & = &  I^T - \frac{2}{v^Tv} \left( vv^T \right)^T  \\
    % & = &  I - \frac{2}{v^Tv} \left( v^T \right)^T v^T \\
    % & = &  I - \frac{2}{v^Tv} vv^T \\
    % & = &  P.
    % \end{eqnarray*}
    % So, $P$ is symmetric.  Now to show that it is orthogonal, we need that $PP^T=I$.  Then
    % \begin{eqnarray*}
    % PP^T & = & P^2  \\
    % & = &  \left( I - \frac{2}{v^Tv} vv^T \right)\left( I - \frac{2}{v^Tv} vv^T \right)  \\
    % & = &  I^2 - \frac{4}{v^Tv} \left( vv^T \right)+ \frac{4}{(v^Tv)^2} \left( vv^Tvv^T \right) \\
    % & = &  I - \frac{4}{v^Tv} vv^T + \frac{4}{v^Tv} \left( vv^T \right) \\
    % & = &  I.
    % \end{eqnarray*}
    % Since $P^T=P$, $P^TP=P^2$ as well, so $P$ is orthogonal.

3. ✍ Let $\mathbf{P}$ be a Householder reflector as in {eq}`hhreflect`.
  
    **(a)** Find a vector $\mathbf{u}$ such that $\mathbf{P}\mathbf{u} = -\mathbf{u}$. ({numref}`fig-hhreflect` may be of help.)

    **(b)** What algebraic condition is necessary and sufficient for a vector $\mathbf{w}$ to satisfy $\mathbf{P}\mathbf{w}=\mathbf{w}$? In $n$ dimensions, how many such linearly independent vectors are there?
  
    (problem-hhsign)=
4. ✍ Under certain circumstances, computing the vector $\mathbf{v}$ in {eq}`hhvector` could lead to subtractive cancellation, which is why line 13 of [`qrfact`](function-qrfact) reads as it does. Devise an example that causes subtractive cancellation if {eq}`hhvector` is used.

5. ✍ Suppose QR factorization is used to compute the solution of a *square* linear system, $\mathbf{A}\mathbf{x}=\mathbf{b}$; i.~e., let $m=n$.

    **(a)** Find an asymptotic flop count for this procedure, and compare to the LU factorization algorithm.

    **(b)** Show that $\kappa_2(\mathbf{A})=\kappa_2(\mathbf{R})$.
  
    % (a) The count for the QR factorization has the largest terms $2mN^2-n^3/3$ as given in
    % the text.  For $m=n$, this reduces to $5n^3/3$ flops.  LU factorization takes $2n^3/3$ flops as stated on
    % p.~55. (Note that the steps for $Q^T b$ and solving the system $Rx=Q^Tb$ are both $O(n^2)$ so that they
    % don't contribute to the asymptotic cost of the system solve by this method.)  So, $QR$ takes more than double the operations for square matrices.
    % (b) Use the definition of the norm.  Note that this approach requires that $A$ be square and invertible.
    % \begin{eqnarray*}
    % ||A||_2 & = &  \max_{||x||_2=1} ||Ax||_2 \\
    % & = &  \max_{||x||_2=1} ||QRx||_2 \\
    % & = &  \max_{||x||_2=1} ||Q(Rx)||_2 \\
    % & = &  \max_{||x||_2=1} ||Rx||_2 \\
    % & = &  ||R||_2. 
    % \end{eqnarray*}
    % Now we're half done.  We also need the same thing for $||A^{-1}||_2$, and now,
    % \begin{eqnarray*}
    % ||A^{-1}||_2 & = &  \max_{||x||_2=1} ||A^{-1}x||_2 \\
    % & = &  \max_{||x||_2=1} ||R^{-1}Q^Tx||_2 \\
    % & = &  \max_{||Qy||_2=1} ||R^{-1}Q^TQy||_2 \\
    % & = &  \max_{||y||_2=1} ||R^{-1}y||_2 \\
    % & = &  ||R^{-1}||_2.
    % \end{eqnarray*}
    % We now have $\kappa_2(A) = ||A||_2||A^{-1}||_2 = ||R||_2||R^{-1}||_2 = \kappa_2(R)$.
    %

6. ✍ Show that $\kappa_2(\mathbf{A})=\kappa_2(\mathbf{R})$ when $\mathbf{A}$ is not square.  (Hint: You can't take an inverse of $\mathbf{A}$ or $\mathbf{R}$.)

    % If we define $\kappa_2(A) = ||A||\, ||A^+||$,
    % then $\kappa_2(R) = ||R||\, ||R^+||$.  Now
    % \[
    % ||A||_2 = ||QR||_2 = ||R||_2
    % \]
    % from our previous homework assignment.  Then
    % \begin{eqnarray*}
    % ||A^+||_2 & = & ||(A^T A)^{-1} A^T||_2 \\
    %          & = & ||(R^T Q^T Q R)^{-1}R^T Q^T||_2 \\
    %          & = & ||(R^T R)^{-1}R^T Q^T||_2 \\
    %          & = & ||(R^T R)^{-1}R^T||_2 \\
    %          & = & ||R^+||_2.
    % \end{eqnarray*}
    % Here we used the fact that an orthogonal matrix does not change the norm, but we never expanded the
    % product inside the norm.  Then
    % \[
    % \kappa_2(A) = ||A||_2 ||A^+||_2 = ||R||_2 ||R^+||_2 = \kappa_2(R).
    % \]
    % Note that we never inverted $A$, $R$ or $Q$ because they are rectangular.

7. ⌨ Modify [`qrfact`](function-qrfact) so that the loops in lines 18--23 are removed. In other words, update the matrices `A` and `Q` with a single statement for each.  Test your code against the original code on an example to verify it.

8. Another algorithmic technique for orthogonally introducing zeros into a matrix is the   *Givens rotation*. Given a 2-vector $[\alpha\, \beta]$, it defines an angle $\theta$ such that
  
    ```{math}
    \begin{bmatrix}
      \cos(\theta) & \sin(\theta) \\ -\sin(\theta) & \cos(\theta)
    \end{bmatrix}
    \begin{bmatrix}
      \alpha \\ \beta
    \end{bmatrix}
    =
    \begin{bmatrix}
      \alpha^2 + \beta^2 \\ 0
    \end{bmatrix}.
    ```

    **(a)** ✍ Given $\alpha$ and $\beta$, show how to compute $\theta$.

    **(b)** ⌨ Given the vector $\mathbf{z}=[1\;2\;3\;4\;5]^T$, use Julia to find a sequence of Givens rotations that transforms $\mathbf{z}$ into the vector $\| \mathbf{z} \|\mathbf{e}_1$. (Hint: You can operate only on pairs of elements at a time, introducing a zero at the lower of the two positions.)
  
