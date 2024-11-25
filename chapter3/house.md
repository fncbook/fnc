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

(section-leastsq-house)=
# Computing QR factorizations

<!--
 It is possible to compute a thin QR factorization using the outer product formula {eq}`matrixouter`, producing the columns and rows of the factors one by one. If $\mathbf{A}=\hat{\mathbf{Q}}\hat{\mathbf{R}}$, then 

 $$
 \hat{\mathbf{Q}}^T\mathbf{A}=\hat{\mathbf{Q}}^T\hat{\mathbf{Q}}\hat{\mathbf{R}} = \hat{\mathbf{R}},
 $$

 by {numref}`Theorem {number} <theorem-qr-ONC>`. Hence if $\mathbf{q}_j$ is the $j$th column of $\hat{\mathbf{Q}}$, then $\mathbf{q}_j\mathbf{A}$ is the $j$th row of $\hat{\mathbf{R}}$.

 (demo-house-qrouter)=
 ```{prf:example}
 ```

 



Since the $\mathbf{q}_j$ have to be unit vectors, we normalize the first column of $\mathbf{A}$ to get $\mathbf{q}_1$.
 
```{code-cell}
A₁ = rand(1:9,6,4)
Q̂,R̂ = zeros(6,4),zeros(4,4)
Q̂[:,1] = normalize(A[:,1]);
```

This leaves us ready to compute the first row of $\hat{\mathbf{R}}$, then subtract off the outer product $\mathbf{q}_1\mathbf{r}_1^T$. 

```{code-cell}
R̂[1,:] = Q̂[:,1]'*A₁
A₂ = A₁ - Q̂[:,1]*R̂[1,:]'
```

This can be repeated.

ETC!




-->

```{index} orthogonal matrix
```

It is possible to compute a thin QR factorization using the outer product formula {eq}`matrixouter`, as we did with LU. However, to stably compute the factorization, a better strategy is to introduce zeros into the lower triangle, one column at a time, using orthogonal matrices. Thanks to {numref}`Theorem %s <theorem-qr-orthogmatrix>`, the product of orthogonal matrices will also be orthogonal.

## Householder reflections

```{index} ! Householder reflector, orthogonal matrix
```

We will use a particular type of orthogonal matrix.

::::{prf:definition} Householder reflector
A **Householder reflector** is a matrix of the form

```{math}
:label: hhreflect
  \mathbf{P} = \mathbf{I} - 2 \mathbf{v} \mathbf{v}^T,
```

where $\mathbf{v}$ is any unit vector (in the 2-norm).
::::

In [Exercise 2](problem-house-reflector) you are asked to show that such a $\mathbf{P}$ is necessarily orthogonal. Note that for any vector $\mathbf{x}$ of appropriate dimension,

```{math}
:label: hhapply
  \mathbf{P}\mathbf{x} = \mathbf{x} - 2 \mathbf{v} (\mathbf{v}^T\mathbf{x}).
```

The reason $\mathbf{P}$ is called a reflector is sketched in {numref}`fig-hhreflect`.

```{figure} figures/hhreflect.svg
:name: fig-hhreflect
A Householder reflector. Because $\mathbf{v}$ is a unit vector, $\mathbf{v}^T\mathbf{x}$ is the component of $\mathbf{x}$ in the direction of $\mathbf{v}$. Hence subtracting $(\mathbf{v}^T\mathbf{x})\mathbf{v}$ projects $\mathbf{x}$ into a hyperplane orthogonal to $\mathbf{v}$. By subtracting off twice as much, we get the reflection of $\mathbf{x}$ through the hyperplane instead.
```

Given a vector $\mathbf{z}$, we can choose $\mathbf{v}$ so that $\mathbf{P}$ reflects $\mathbf{z}$ onto the $x_1$-axis—i.e., so that $\mathbf{P}\mathbf{z}$ is nonzero only in the first element. Because orthogonal matrices preserve the 2-norm, we must have

```{math}
:label: hhgoal
\mathbf{P}\mathbf{z} =
\begin{bmatrix}
\pm \| \mathbf{z} \|\\0 \\ \vdots \\ 0
\end{bmatrix} = \pm \| \mathbf{z} \| \mathbf{e}_1.
```

(Recall that $\mathbf{e}_k$ is the $k$th column of the identity matrix.) We choose the positive sign above for our discussion, but see {numref}`Function {number} <function-qrfact>` and [Exercise 4](problem-house-sign) for important computational details. Let

```{math}
:label: hhvector
  \mathbf{w} = \| \mathbf{z} \| \mathbf{e}_1-\mathbf{z}, \quad \mathbf{v} = \frac{\mathbf{w}}{\|\mathbf{w}\|}.
```

If it turns out that $\mathbf{w}=\boldsymbol{0}$, then $\mathbf{z}$ is already in the target form and we can take $\mathbf{P}=\mathbf{I}$.  Otherwise, we have the following.

(theorem-hhreflect)=
````{prf:theorem} Householder reflector
Let $\mathbf{v}$ be defined by {eq}`hhvector` and let $\mathbf{P}$ be given by {eq}`hhreflect`. Then $\mathbf{P}$ is symmetric and orthogonal, and $\mathbf{P}\mathbf{z}=\| \mathbf{z} \|\mathbf{e}_1$.
````

````{prf:proof}
The proofs of symmetry and orthogonality are left to [Exercise 2](problem-house-reflector). For the last fact, we use {eq}`hhapply` to compute
  
```{math}
\mathbf{P}\mathbf{z} = \mathbf{z} - 2 \frac{\mathbf{w}^T \mathbf{z}}{\mathbf{w}^T\mathbf{w}} \mathbf{w}.
```

Since $\mathbf{e}_1^T\mathbf{z}=z_1$,
  
```{math}
\begin{split}
    \mathbf{w}^T\mathbf{w} &= \| \mathbf{z} \|^2 - 2 \| \mathbf{z} \| z_1 + \mathbf{z}^T\mathbf{z}
    = 2\| \mathbf{z} \|(\| \mathbf{z} \|-z_1),\\
    \mathbf{w}^T\mathbf{z} &= \| \mathbf{z} \|z_1 - \mathbf{z}^T\mathbf{z} = -\| \mathbf{z} \|\bigl(\| \mathbf{z} \|-z_1\bigr),
\end{split}
```

leading finally to
  
```{math}
\mathbf{P}\mathbf{z} = \mathbf{z} - 2\cdot
\frac{-\| \mathbf{z} \| \bigl(\| \mathbf{z} \|-z_1\bigr)}{2\| \mathbf{z} \| \bigl(\| \mathbf{z} \|-z_1\bigr)} \mathbf{w}
= \mathbf{z} + \mathbf{w} = \| \mathbf{z} \|\mathbf{e}_1.
```
````

## Factorization algorithm

```{index} matrix factorization; QR
```

The QR factorization is computed by using successive Householder reflections to introduce zeros in one column at a time. We first show the process for a small numerical example in {numref}`Demo %s <demo-house-qr>`.

(demo-house-qr)=
```{prf:example}
```





::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item}
:columns: 7


We will use Householder reflections to produce a QR factorization of a random matrix.


:::
:::{grid-item-card}
:columns: 5


The `rand` function can select randomly from within the interval $[0,1]$, or from a vector or range that you specify.

:::
::::

```{code-cell}
A = rand(float(1:9),6,4)
m,n = size(A)
```

```{index} Julia; normalize, ! Julia; I
```

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item}
:columns: 7


Our first step is to introduce zeros below the diagonal in column 1 by using {eq}`hhvector` and {eq}`hhreflect`. 


:::
:::{grid-item-card}
:columns: 5


`I` can stand for an identity matrix of any size, inferred from the context when needed.

:::
::::

```{code-cell}
z = A[:,1];
v = normalize(z - norm(z)*[1;zeros(m-1)])
P₁ = I - 2v*v'   # reflector
```

We check that this reflector introduces zeros as it should:

```{code-cell}
P₁*z
```

Now we replace $\mathbf{A}$ by $\mathbf{P}\mathbf{A}$.

```{code-cell}
A = P₁*A
```

We are set to put zeros into column 2. We must not use row 1 in any way, lest it destroy the zeros we just introduced. So we leave it out of the next reflector.

```{code-cell}
z = A[2:m,2]
v = normalize(z - norm(z)*[1;zeros(m-2)])
P₂ = I - 2v*v'
```

We now apply this reflector to rows 2 and below only.

```{code-cell}
A[2:m,:] = P₂*A[2:m,:]
A
```

We need to iterate the process for the last two columns.

```{code-cell}
for j in 3:n
    z = A[j:m,j]
    v = normalize(z - norm(z)*[1;zeros(m-j)])
    P = I - 2v*v'
    A[j:m,:] = P*A[j:m,:]
end
```

We have now reduced the original to an upper triangular matrix using four orthogonal Householder reflections:

```{code-cell}
R = triu(A)
```




You may be wondering what happened to $\mathbf{Q}$ in {numref}`Demo %s <demo-house-qr>`. Each Householder reflector is orthogonal but not full-size. We have to pad it out to represent algebraically the fact that a block of the first rows is left alone. Given a reflector $\mathbf{P}_k$ that is of square size $m-k+1$, we define

```{math}
\mathbf{Q}_k =
\begin{bmatrix}
\mathbf{I}_{k-1} & \boldsymbol{0} \\ \boldsymbol{0} & \mathbf{P}_k
\end{bmatrix}.
```

It is easy to show that $\mathbf{Q}_k$ is also orthogonal. Then the algorithm produces

```{math}
  \mathbf{Q}_n \mathbf{Q}_{n-1}\cdots \mathbf{Q}_1 \mathbf{A} = \mathbf{R}.
```

But $\mathbf{Q}_n \mathbf{Q}_{n-1}\cdots \mathbf{Q}_1$ is orthogonal too, and we multiply on the left by its transpose to get $\mathbf{A}=\mathbf{Q}\mathbf{R}$, where $\mathbf{Q} =  (\mathbf{Q}_n \mathbf{Q}_{n-1}\cdots \mathbf{Q}_1)^T$. We don't even need to form these matrices explicitly. Writing

```{math}
:label: hhbuildQ
\mathbf{Q}^T = \mathbf{Q}_n \mathbf{Q}_{n-1}\cdots \mathbf{Q}_1 = \mathbf{Q}_n \Bigl( \mathbf{Q}_{n-1}\bigl(\cdots (\mathbf{Q}_1\mathbf{I})\cdots\bigr)\Bigr),
```

we can build $\mathbf{Q}^T$ iteratively by starting with the identity and doing the same row operations as on $\mathbf{A}$. That process uses much less memory than building the $\mathbf{Q}_k$ matrices explicitly.

The algorithm we have described is encapsulated in {numref}`Function {number} <function-qrfact>`. There is one more refinement in it, however. As indicated by {eq}`hhapply`, the application of a reflector $\mathbf{P}$ to a vector does not require the formation of the matrix $\mathbf{P}$ explicitly.

(function-qrfact)=

```{prf:function} qrfact
**QR factorization by Householder reflections**
```{code-block} julia
:lineno-start: 1
"""
    qrfact(A)

QR factorization by Householder reflections. Returns Q and R.
"""
function qrfact(A)
    m,n = size(A)
    Qt = diagm(ones(m))
    R = float(copy(A))
    for k in 1:n
        z = R[k:m,k]
        w = [ -sign(z[1])*norm(z) - z[1]; -z[2:end] ]
        nrmw = norm(w)
        if nrmw < eps() continue; end    # skip this iteration
        v = w / nrmw;
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

## Q-less QR and least squares

In {numref}`Demo {number} <demo-qr-qrfact>` it was seen that the $\mathbf{Q}$ output of Julia's `qr` function is not a standard matrix. The reason is that Equation {eq}`lsqr` shows that in order to solve the linear least-squares problem, all we need from $\mathbf{Q}$ is the computation of $\hat{\mathbf{Q}}^T\mathbf{b}$. Referring again to {eq}`hhbuildQ` and {eq}`hhapply`, the special structure of the reflectors is such that for this computation, we only need to apply code similar to lines 18 and 21 of {numref}`Function {number} <function-qrfact>` for each of the Householder vectors $\mathbf{v}$ that is constructed. 

This observation leads to the idea of the *Q-less QR factorization*, in which the full or thin $\mathbf{Q}$ is never computed explicitly. This is the variant used by Julia's `qr`. The returned value `Q` used within {numref}`Function {number} <function-lsqrfact>` is of a special type that allows Julia to perform `Q'*b` efficiently for the least-squares solution. 

In [Exercise 8](problem-house-flops) you are asked to derive the following result about the Q-less factorization.

(theorem-house-flops)=
:::{prf:theorem} 
Q-less QR factorization by Householder reflections takes $\sim(2mn^2-\frac{2}{3}n^3)$ flops.
:::

The flop count quoted in {numref}`Theorem {number} <theorem-house-flops>` dominates the running time for least-squares solution via QR. Compared to the count from {numref}`Theorem {number} <theorem-normaleqns-flops>` for solution by the normal equations, the flops are essentially identical when $m=n$, but the QR solution is about twice the cost when $m\gg n$. The redeeming quality of the QR route is better stability, which we do not discuss here.

## Exercises

1. ⌨ Find a Householder reflector $\mathbf{P}$ such that
  
    ```{math}
    \mathbf{P}
    \begin{bmatrix}
      2 \\ 9 \\ -6
    \end{bmatrix} =
    \begin{bmatrix}
      11\\0\\0
    \end{bmatrix}.
    ```
    (problem-house-reflector)=
2. ✍ Prove the unfinished items in {numref}`Theorem %s <theorem-hhreflect>`, namely that a Householder reflector $\mathbf{P}$ is symmetric and orthogonal.

3. ✍ Let $\mathbf{P}$ be a Householder reflector as in {eq}`hhreflect`.
  
    **(a)** Find a vector $\mathbf{u}$ such that $\mathbf{P}\mathbf{u} = -\mathbf{u}$. ({numref}`fig-hhreflect` may be of help.)

    **(b)** What algebraic condition is necessary and sufficient for a vector $\mathbf{x}$ to satisfy $\mathbf{P}\mathbf{x}=\mathbf{x}$? In $n$ dimensions, how many such linearly independent vectors are there?
  
    (problem-house-sign)=
4. ✍ Under certain circumstances, computing the vector $\mathbf{v}$ in {eq}`hhvector` could lead to subtractive cancellation, which is why line 12 of {numref}`Function {number} <function-qrfact>` reads as it does. Devise an example that causes subtractive cancellation if {eq}`hhvector` is used.

5. ✍ Suppose QR factorization is used to compute the solution of a *square* linear system, $\mathbf{A}\mathbf{x}=\mathbf{b}$, i.e., let $m=n$.

    **(a)** Find an asymptotic flop count for this procedure, and compare it to the LU factorization algorithm.

    **(b)** Show that $\kappa_2(\mathbf{A}) = \kappa_2(\mathbf{R})$.
  
6. ✍ Prove that $\kappa_2(\mathbf{A})=\kappa_2(\mathbf{R})$ when $\mathbf{A}$ is not square.  (Be careful! You can't take an inverse of $\mathbf{A}$ or $\mathbf{R}$.)

7. Another algorithmic technique for orthogonally introducing zeros into a matrix is the   *Givens rotation*. Given a 2-vector $[\alpha,\, \beta]$, it defines an angle $\theta$ such that
  
    ```{math}
    \begin{bmatrix}
      \cos(\theta) & \sin(\theta) \\ -\sin(\theta) & \cos(\theta)
    \end{bmatrix}
    \begin{bmatrix}
      \alpha \\ \beta
    \end{bmatrix}
    =
    \begin{bmatrix}
      \sqrt{\alpha^2 + \beta^2} \\ 0
    \end{bmatrix}.
    ```

    **(a)** ✍ Given $\alpha$ and $\beta$, show how to compute $\cos \theta$ and $\sin \theta$ without evaluating any trig functions.

    **(b)** ⌨ Given the vector $\mathbf{z}=[1\;2\;3\;4\;5]^T$, use Julia to find a sequence of Givens rotations that transforms $\mathbf{z}$ into the vector $\| \mathbf{z} \|\mathbf{e}_1$. (Hint: You can operate only on pairs of elements at a time, introducing a zero at the lower of the two positions.)
  
    (problem-house-flops)=
8. ✍ Derive the result of {numref}`Theorem {number} <theorem-house-flops>` by analyzing {numref}`Function {number} <function-qrfact>` without lines 20–22. 

    (problem-house-speed)=
9. ✍ Suppose $m=Kn$ for constant $K \ge 1$ as both $m$ and $n$ go to infinity. Show that the flop counts from {numref}`Theorem {number} <theorem-house-flops>`  and {numref}`Theorem {number} <theorem-normaleqns-flops>` have a ratio of 1 when $K=1$ and approaches 2 as $K\to \infty$. 
