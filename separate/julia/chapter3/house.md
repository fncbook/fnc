---
numbering:
  enumerator: 3.4.%s
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.12
---
```{code-cell}
:tags: [remove-cell]
import Pkg; Pkg.activate("/Users/driscoll/Documents/GitHub/fnc")

using FNCFunctions

using Plots
default(
    titlefont=(11,"Helvetica"),
    guidefont=(11,"Helvetica"),
    linewidth = 2,
    markersize = 3,
    msa = 0,
    size=(500,320),
    label="",
    html_output_format = "svg"
)

using PrettyTables, LaTeXStrings, Printf
using LinearAlgebra

@ptconf backend = Val(:html) tf = tf_html_simple
```

(section-leastsq-house)=

# Computing QR factorizations

```{index} orthogonal matrix
```

It is possible to compute a thin QR factorization using the outer product formula {eq}`matrixouter`, as we did with LU. However, to stably compute the factorization, a better strategy is to introduce zeros into the lower triangle, one column at a time, using orthogonal matrices. Thanks to @theorem-qr-orthogmatrix, the product of orthogonal matrices will also be orthogonal.

## Householder reflections

We will use a particular type of orthogonal matrix.

```{index} ! Householder reflector, orthogonal matrix
```

::::{prf:definition} Householder reflector
:label: definition-hhreflect
A **Householder reflector** is a matrix of the form

```{math}
:label: hhreflect
  \mathbf{P} = \mathbf{I} - 2 \mathbf{v} \mathbf{v}^T,
```

where $\mathbf{v}$ is any unit vector (in the 2-norm).
::::

In @problem-house-reflector you are asked to show that such a $\mathbf{P}$ is necessarily orthogonal. Note that for any vector $\mathbf{x}$ of appropriate dimension,

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

(Recall that $\mathbf{e}_k$ is the $k$th column of the identity matrix.) We choose the positive sign above for our discussion, but see {numref}`Function {number} <function-qrfact>` and @problem-house-sign for important computational details. Let

```{math}
:label: hhvector
  \mathbf{w} = \| \mathbf{z} \| \mathbf{e}_1-\mathbf{z}, \quad \mathbf{v} = \frac{\mathbf{w}}{\|\mathbf{w}\|}.
```

If it turns out that $\mathbf{w}=\boldsymbol{0}$, then $\mathbf{z}$ is already in the target form and we can take $\mathbf{P}=\mathbf{I}$.  Otherwise, we have the following.

````{prf:theorem} Householder reflector
:label: theorem-hhreflect
Let $\mathbf{v}$ be defined by {eq}`hhvector` and let $\mathbf{P}$ be given by {eq}`hhreflect`. Then $\mathbf{P}$ is symmetric and orthogonal, and $\mathbf{P}\mathbf{z}=\| \mathbf{z} \|\mathbf{e}_1$.
````

````{prf:proof}
:enumerated: false

The proofs of symmetry and orthogonality are left to @problem-house-reflector. For the last fact, we use {eq}`hhapply` to compute
  
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

The QR factorization is computed by using successive Householder reflections to introduce zeros in one column at a time. We first show the process for a small numerical example in @demo-house-qr.

::::{prf:example} Householder QR factorization
:label: demo-house-qr


We will use Householder reflections to produce a QR factorization of a random matrix.
```{tip}
:class: dropdown
The `rand` function can select randomly from within the interval $[0,1]$, or from a vector or range that you specify.
```

```{code-cell}
A = rand(float(1:9), 6, 4)
m,n = size(A)
```

```{index} Julia; normalize, ! Julia; I
```

Our first step is to introduce zeros below the diagonal in column 1 by using {eq}`hhvector` and {eq}`hhreflect`. 

```{tip}
:class: dropdown
`I` can stand for an identity matrix of any size, inferred from the context when needed.
```

```{code-cell}
z = A[:, 1];
v = normalize(z - norm(z) * [1; zeros(m-1)])
P₁ = I - 2v * v'   # reflector
```

We check that this reflector introduces zeros as it should:

```{code-cell}
P₁ * z
```

Now we replace $\mathbf{A}$ by $\mathbf{P}\mathbf{A}$.

```{code-cell}
A = P₁ * A
```

We are set to put zeros into column 2. We must not use row 1 in any way, lest it destroy the zeros we just introduced. So we leave it out of the next reflector.

```{code-cell}
z = A[2:m, 2]
v = normalize(z - norm(z) * [1; zeros(m-2)])
P₂ = I - 2v * v'
```

We now apply this reflector to rows 2 and below only.

```{code-cell}
A[2:m, :] = P₂ * A[2:m, :]
A
```

We need to iterate the process for the last two columns.

```{code-cell}
for j in 3:n
    z = A[j:m, j]
    v = normalize(z - norm(z) * [1; zeros(m-j)])
    P = I - 2v * v'
    A[j:m, :] = P * A[j:m, :]
end
```

We have now reduced the original to an upper triangular matrix using four orthogonal Householder reflections:

```{code-cell}
R = triu(A)
```

::::

You may be wondering what happened to $\mathbf{Q}$ in @demo-house-qr. Each Householder reflector is orthogonal but not full-size. We have to pad it out to represent algebraically the fact that a block of the first rows is left alone. Given a reflector $\mathbf{P}_k$ that is of square size $m-k+1$, we define

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

``````{prf:algorithm} qrfact
:label: function-qrfact

```{literalinclude} chapter03.jl
:filename: qrfact.jl
:language: julia
:start-after: # begin qrfact
:end-before: # end qrfact
:linenos: true
```
``````

## Q-less QR and least squares

In @demo-qr-qrfact, it was seen that the $\mathbf{Q}$ output of Julia's `qr` function is not a standard matrix. The reason is that Equation {eq}`lsqr` shows that in order to solve the linear least-squares problem, all we need from $\mathbf{Q}$ is the computation of $\hat{\mathbf{Q}}^T\mathbf{b}$. Referring again to {eq}`hhbuildQ` and {eq}`hhapply`, the special structure of the reflectors is such that for this computation, we only need to apply code similar to lines 18 and 21 of {numref}`Function {number} <function-qrfact>` for each of the Householder vectors $\mathbf{v}$ that is constructed.

This observation leads to the idea of the *Q-less QR factorization*, in which the full or thin $\mathbf{Q}$ is never computed explicitly. This is the variant used by Julia's `qr`. The returned value `Q` used within {numref}`Function {number} <function-lsqrfact>` is of a special type that allows Julia to perform `Q'*b` efficiently for the least-squares solution.

In @problem-house-flops you are asked to derive the following result about the Q-less factorization.

:::{prf:theorem}
:label: theorem-house-flops
Q-less QR factorization by Householder reflections takes $\sim(2mn^2-\frac{2}{3}n^3)$ flops.
:::

The flop count quoted in @theorem-house-flops dominates the running time for least-squares solution via QR. Compared to the count from @theorem-normaleqns-flops for solution by the normal equations, the flops are essentially identical when $m=n$, but the QR solution is about twice the cost when $m\gg n$. The redeeming quality of the QR route is better stability, which we do not discuss here.

## Exercises

``````{exercise}
:label: problem-house-find
⌨ Find a Householder reflector $\mathbf{P}$ such that

```{math}
:numbered: false
\mathbf{P}
\begin{bmatrix}
2 \\ 9 \\ -6
\end{bmatrix} =
\begin{bmatrix}
11\\0\\0
\end{bmatrix}.
```
``````

``````{exercise}
:label: problem-house-reflector
✍ Prove the unfinished items in @theorem-hhreflect, namely, that a Householder reflector $\mathbf{P}$ is symmetric and orthogonal.
``````

``````{exercise}
:label: problem-house-negate
✍ Let $\mathbf{P}$ be a Householder reflector as in {eq}`hhreflect`.

**(a)** Find a vector $\mathbf{u}$ such that $\mathbf{P}\mathbf{u} = -\mathbf{u}$. (@fig-hhreflect may be of help.)

**(b)** What algebraic condition is necessary and sufficient for a vector $\mathbf{x}$ to satisfy $\mathbf{P}\mathbf{x}=\mathbf{x}$? In $n$ dimensions, how many such linearly independent vectors are there?

``````

``````{exercise}
:label: problem-house-sign
✍ Under certain circumstances, computing the vector $\mathbf{v}$ in {eq}`hhvector` could lead to subtractive cancellation, which is why `w` is computed in a particular way in @function-qrfact. Devise an example that causes subtractive cancellation if {eq}`hhvector` is used.
``````

``````{exercise}
:label: problem-house-square
✍ Suppose QR factorization is used to compute the solution of a *square* linear system, $\mathbf{A} \mathbf{x}=\mathbf{b}$, i.e., let $m=n$.

**(a)** Find an asymptotic flop count for this procedure, and compare it to the LU factorization algorithm.

**(b)** Show that $\kappa_2(\mathbf{A}) = \kappa_2(\mathbf{R})$.
``````

``````{exercise}
:label: problem-house-cond
✍ Prove that $\kappa_2(\mathbf{A})=\kappa_2(\mathbf{R})$ when $\mathbf{A}$ is not square.  (Be careful! You can't take an inverse of $\mathbf{A}$ or $\mathbf{R}.$)
``````

``````{exercise}
:label: problem-house-givens
Another algorithmic technique for orthogonally introducing zeros into a matrix is the *Givens rotation*. Given a 2-vector $[\alpha,\, \beta],$ it defines an angle $\theta$ such that

```{math}
:numbered: false
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

**(b)** ⌨ Given the vector $\mathbf{z}=[1\;2\;3\;4\;5]^T,$ find a sequence of Givens rotations that transforms $\mathbf{z}$ into the vector $\| \mathbf{z} \|\mathbf{e}_1.$ (Hint: You can operate only on pairs of elements at a time, working from the bottom upwards.)

``````

``````{exercise}
:label: problem-house-flops
✍ Derive the result of @theorem-house-flops by analyzing {numref}`Function {number} <function-qrfact>`, excluding the loop that updates $\mathbf{Q}.$
``````

``````{exercise}
:label: problem-house-speed
✍ Suppose $m=Kn$ for constant $K \ge 1$ as both $m$ and $n$ go to infinity. Show that the flop counts from @theorem-house-flops and @theorem-normaleqns-flops have a ratio of 1 when $K=1$ and approaches 2 as $K\to \infty.$ 
``````
