# Krylov subspaces

```{index} matrix; Krylov
```
```{index} Krylov matrix
```

The power and inverse iterations have a flaw that seems obvious once it is pointed out. Given a seed vector $\mathbf{u}$, they produce a sequence of vectors $\mathbf{u}_1,\mathbf{u}_2,\ldots$ that are scalar multiples of $\mathbf{u},\mathbf{A}\mathbf{u},\mathbf{A}^{2}\mathbf{u},\ldots$, but only the most recent vector is used to produce an eigenvector estimate. It stands to reason that we could do no worse, and perhaps much better, if we searched among all linear combinations of the vectors seen in the past. In other words, we seek a solution in the range (column space) of the $m\times n$ {term}`Krylov matrix`

:::{math}
:label: krylovmatrix
\mathbf{K}_m =
\begin{bmatrix}
  \mathbf{u} & \mathbf{A}\mathbf{u} & \mathbf{A}^{2} \mathbf{u} & \cdots & \mathbf{A}^{m-1} \mathbf{u}
\end{bmatrix}.
:::

```{index} Krylov subspace
```
Such a space is called the $m$th {term}`Krylov subspace` $\mathcal{K}_m$ of $\mathbb{C}^n$.\footnote{A proper pronunciation of "Krylov" is something like "kree-luv," but American English speakers often say "kreye-lahv."} Implicitly we understand that $\mathbf{K}_m$ and $\mathcal{K}_m$ depend on both $\mathbf{A}$ and the initial vector $\mathbf{u}$, but we rarely express the dependence notationally. In general, we expect that the dimension of $\mathcal{K}_m$, which is the rank of $\mathbf{K}_m$, equals $m$, though it may be smaller.

## Properties

As we have seen with the power iteration, part of the appeal of the Krylov matrix is that it can be generated in a way that fully exploits the sparsity of $\mathbf{A}$, simply through repeated matrix-vector multiplication. Furthermore, we have some important mathematical properties.

::::{prf:theorem}
:label: theorem-krylovmult
Suppose $\mathbf{A}$ is $n\times n$, $0<m<n$, and a vector $\mathbf{u}$ is used to generate Krylov subspaces. If $\mathbf{x}\in\mathcal{K}_m$, then the following hold:

1. $\mathbf{x} = \mathbf{K}_m \mathbf{z}$ for some $\mathbf{z}\in\mathbb{C}^m$.
2. $\mathbf{x} \in \mathcal{K}_{m+1}$.
3. $\mathbf{A}\mathbf{x} \in \mathcal{K}_{m+1}$.
::::

::::{prf:proof}
If $\mathbf{x}\in\mathcal{K}_m$, then for some coefficients $c_1,\ldots,c_m$,

:::{math}
\mathbf{x} = c_1 \mathbf{u} + c_2 \mathbf{A} \mathbf{u} + \cdots + c_m \mathbf{A}^{m-1} \mathbf{u}.
:::

Thus let $\mathbf{z}= \begin{bmatrix} c_1 & \cdots & c_m \end{bmatrix}^T$. Also $\mathbf{x}\in\mathcal{K}_{m+1}$, as we can add zero times $\mathbf{A}^{m}\mathbf{u}$ to the sum. Finally,
  
:::{math}
\mathbf{A}\mathbf{x} = c_1 \mathbf{A} \mathbf{u} + c_2 \mathbf{A}^{2} \mathbf{u} + \cdots + c_m \mathbf{A}^{m} \mathbf{u} \in \mathcal{K}_{m+1}.
:::
::::

## Reducing dimension

```{index} dimension reduction
```
The problems $\mathbf{A}\mathbf{x}=\mathbf{b}$ and $\mathbf{A}\mathbf{x}=\lambda\mathbf{x}$ are statements about a very high-dimensional space $\mathbb{C}^n$. One way to approximate them is to replace the full $n$-dimensional space with a much lower-dimensional $\mathcal{K}_m$ for $m\ll n$. This is the essence of the Krylov subspace approach.

For instance, we can interpret $\mathbf{A}\mathbf{x}_m\approx \mathbf{b}$ in the sense of linear least squares—that is, using \lemref{krylovmult} to let $\mathbf{x}=\mathbf{K}_m\mathbf{z}$,

:::{math}
:label: gmresdef
\min_{\mathbf{x}\in\\\mathcal{K}_m} \|  \mathbf{A}\mathbf{x}-\mathbf{b} \|
= \min_{\mathbf{z}\in\mathbb{C}^m} \| \mathbf{A}(\mathbf{K}_m\mathbf{z})-\mathbf{b} \|
= \min_{\mathbf{z}\in\mathbb{C}^m} \| (\mathbf{A}\mathbf{K}_m)\mathbf{z})-\mathbf{b} \|.
:::

The natural seed vector for $\mathcal{K}_m$ in this case is the vector $\mathbf{b}$. In the next example we try to implement {eq}`gmresdef`. We do take one precaution: because the vectors $\mathbf{A}^{k}\mathbf{b}$ may become very large or small in norm, we normalize after each multiplication by $\mathbf{A}$, just as we did in the power iteration.

(demo-subspace-unstable)=
::::{prf:example} Julia demo
:class: demo
:label: demos-subspace-unstable
{doc}`demos/subspace-unstable`
::::

## The Arnoldi iteration

The breakdown of convergence in {prf:ref}`demos-subspace-unstable` is due to a critical numerical defect in our approach: the columns of the Krylov matrix {eq}`krylovmatrix` increasingly become parallel to one another (and the dominant eigenvector), as {eq}`poweriterconverge` predicts. As we saw [previously](../leastsq/qr.md), near-parallel vectors create the potential for numerical cancellation. This manifests as a large condition number for $\mathbf{K}_m$ as $m$ grows, eventually creating excessive error when solving the least-squares system.

The polar opposite of an ill-conditioned basis for $\mathcal{K}_m$ is an orthonormal one. Suppose we had a thin QR factorization of $\mathbf{K}_m$:

\begin{align*}
  \mathbf{K}_m  = \mathbf{Q}_m \mathbf{R}_m
  & =
  \begin{bmatrix}
    \mathbf{q}_1& \mathbf{q}_2 & \cdots & \mathbf{q}_m
  \end{bmatrix}
  \begin{bmatrix}
    R_{11} & R_{12} & \cdots & R_{1m} \\
    0 & R_{22} & \cdots & R_{2m} \\
    \vdots & & \ddots & \\
    0 & 0 & \cdots & R_{mm}
  \end{bmatrix}.
\end{align*}

Then the vectors $\mathbf{q}_1$, \ldots, $\mathbf{q}_m$ are the orthonormal basis we seek for $\mathcal{K}_m$. By {prf:ref}`theorem-krylovmult`, we know that $\mathbf{A}\mathbf{q}_m \in \mathcal{K}_{m+1}$, and therefore

:::{math}
:label: arnoldivec
\mathbf{A} \mathbf{q}_m = H_{1m} \, \mathbf{q}_1 + H_{2m} \, \mathbf{q}_2 + \cdots + H_{m+1,m}\,\mathbf{q}_{m+1},
:::

for some choice of the $H_{ij}$. Note that by using orthonormality, we have

:::{math}
:label: arnoldiip
\mathbf{q}_i^* (\mathbf{A}\mathbf{q}_m) = H_{im},\qquad i=1,\ldots,m.
:::

Since we started by assuming that we know $\mathbf{q}_1,\ldots,\mathbf{q}_m$, the only unknowns in {eq}`arnoldivec` are $H_{m+1,m}$ and $\mathbf{q}_{m+1}$. But they appear only as a product, and we know that $\mathbf{q}_{m+1}$ is a *unit* vector, so they are uniquely defined (up to sign) by the other terms in the equation.

We can now proceed iteratively. If $\mathbf{u}$ is the Krylov seed vector, then

1. Let $\mathbf{q}_1=\mathbf{u}/\|\mathbf{u}\|$.
1. For $m=1,2,\ldots$
    1. Use {eq}`arnoldiip` to find $H_{im}$ for $i=1,\ldots,m$.
    1. Let

    :::{math}
    :label: arnoldigs
    \mathbf{v} = (\mathbf{A} \mathbf{q}_m) - H_{1m} \,\mathbf{q}_1 - H_{2m}\, \mathbf{q}_2 - \cdots - H_{mm}\, \mathbf{q}_m.
    :::

    1. Let $H_{m+1,m}=\|\mathbf{v}\|$.
    1. Let $\mathbf{q}_{m+1}=\mathbf{v}/H_{m+1,m}$.

```{index} Arnoldi iteration
```
We have just described the {term}`Arnoldi iteration`. The Arnoldi iteration finds  orthonormal bases for a nested sequence of Krylov subspaces.

::::{prf:example} Julia demo
:class: demo
:label: demos-subspace-arnoldi
{doc}`demos/subspace-arnoldi`
::::

## Key identity

Up to now we have focused only on finding the orthonormal basis that lies in the columns of $\mathbf{Q}_m$. But the $H_{ij}$ values found during the iteration are also important. Taking $j=1,2,\ldots,m$ in {eq}`arnoldivec` leads to

:::{math}
:label: arnoldimat
\begin{split}
  \mathbf{A}\mathbf{Q}_m &= \begin{bmatrix}
    \mathbf{A}\mathbf{q}_1 & \cdots \mathbf{A}\mathbf{q}_m
  \end{bmatrix}\\
  & = \begin{bmatrix}
    \mathbf{q}_1 & \mathbf{q}_2 & \cdots & \mathbf{q}_{m+1}
  \end{bmatrix} \begin{bmatrix}
    H_{11} & H_{12} & \cdots & H_{1m} \\
    H_{21} & H_{22} & \cdots & H_{2m} \\
    & H_{32} & \ddots & \vdots \\
    & & \ddots & H_{mm} \\
    & & & H_{m+1,m}
\end{bmatrix} = \mathbf{Q}_{m+1} \mathbf{H}_m,
\end{split}
:::

```{index} matrix; upper Hessenberg
```
where the matrix $\mathbf{H}_m$ has an {term}`upper Hessenberg` structure. Equation {eq}`arnoldimat` is a fundamental identity of Krylov subspace methods.

## Implementation

(function-arnoldi)=
````{proof:function} arnoldi
**Arnoldi iteration for Krylov subspaces**

```{code-block} julia
:lineno-start: 1
"""
arnoldi(A,u,m)

Perform the Arnoldi iteration for `A` starting with vector `u`, out
to the Krylov subspace of degree `m`. Return the orthonormal basis
(`m`+1 columns) and the upper Hessenberg `H` of size `m`+1 by `m`.
"""
function arnoldi(A,u,m)
    n = length(u)
    Q = zeros(n,m+1)
    H = zeros(m+1,m)
    Q[:,1] = u/norm(u)
    for j = 1:m
      # Find the new direction that extends the Krylov subspace.
      v = A*Q[:,j]
      # Remove the projections onto the previous vectors.
      for i = 1:j
        H[i,j] = dot(Q[:,i],v)
        v -= H[i,j]*Q[:,i]
      end
      # Normalize and store the new basis vector.
      H[j+1,j] = norm(v)
      Q[:,j+1] = v/H[j+1,j]
    end

    return Q,H
end
"""
```
````

An implementation of the Arnoldi iteration is given in {numref}`Function {number}<function-arnoldi>`. A careful inspection shows that the loop at line~18 does not exactly implement {eq}`arnoldiip` and {eq}`arnoldigs`. The reason is numerical stability. Though the described and implemented versions are mathematically equivalent in exact arithmetic (see {ref}`prob-arnoldi-modifiedgs`), the approach in {ref}`fun-arnoldi` is much more stable to roundoff.

In the next section we revisit the idea of approximately solving $\mathbf{A}\mathbf{x}=\mathbf{b}$ over a Krylov subspace $\mathcal{K}_m$, using the ONC matrix $\mathbf{Q}_m$ in place of $\mathbf{K}_m$. A [related idea](`problem-krylov-arnoldieig`) is used to approximate the eigenvalue problem for $\mathbf{A}$, which is the approach that underlies `eigs` for sparse matrices.

## Exercises

(problem-krylovpermute)=
1. ✍ Consider the matrix $\mathbf{A}=\displaystyle \begin{bmatrix}
    0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \\ 1 & 0 & 0 & 0
    \end{bmatrix}.$

    **(a)** Find the Krylov matrix $\mathbf{K}_3$ for the seed vector $\mathbf{u}=\mathbf{e}_1$. 

    **(b)** Find $\mathbf{K}_3$ for the seed vector $\mathbf{u}=\begin{bmatrix}1 & 1 & 1 & 1\end{bmatrix}^T.$  

    ::::{only} solutions
        \begin{enumerate}
        **()** The Krylov vectors are $b=e_1$, $Ab=e_4$, $A^2b=e_3$,
          $A^3b=e_2$. 
        **()** Application of $\mathbf{A}$ does not change the vector. 
        \end{enumerate}
    ::::

2. ⌨ For each matrix, make a table of the 2-norm condition numbers $\kappa(\mathbf{K}_m)$ for $m=1,\ldots,10$. Use a vector of all ones as the Krylov seed.

    **(a)** Matrix from {prf:ref}`demos-subspace-unstable`
    
    **(b)** $\begin{bmatrix}
      -2 & 1 & & &  \\
      1 & -2 & 1 & &  \\
      & \ddots & \ddots & \ddots & \\
      & & 1 & -2 & 1 \\
       & & & 1 & -2
     \end{bmatrix} \: (100\times 100)$

   **(c)** $\begin{bmatrix}
      -2 & 1 & & & 1 \\
      1 & -2 & 1 & &  \\
      & \ddots & \ddots & \ddots & \\
      & & 1 & -2 & 1 \\
      1 & & & 1 & -2
    \end{bmatrix} \:(200\times 200) $

    ::::{only} solutions
    lambda = 10 + (1:100);
    A = diag(lambda) + triu(rand(100),1);
    u = ones(100,1);
    K = [];
    for m = 1:10
        K = [ K A^(m-1)*u ];
        kappa(m) = cond(K);
    end
    table(kappa(:))    
    ::::

    (problem-matrixpolykrylov)=
3. ✍ (See also an [earlier exercise](problem-matrixpolyeig).) If $p(z)=c_0+c_1z+\cdots+c_kz^k$ is a polynomial, then its value for a matrix argument is defined as
  
    $$
    p(\mathbf{A})= c_0\mathbf{I} + c_1 \mathbf{A} + \cdots +c_k \mathbf{A}^k.
    $$

    Show that if $\mathbf{x}\in\mathcal{K}_m$, then $\mathbf{x}=p(\mathbf{A})\mathbf{u}$ for a polynomial $p$ of degree at most $m-1$. 

4. ✍ Compute the asymptotic flop requirements for {numref}`Function {number}<function-arnoldi>`. Assume that due to sparsity, a matrix-vector multiplication $\mathbf{A}\mathbf{u}$ requires only $c n$ flops for a constant $c$, rather than the usual $O(n^2)$. 

5. ⌨ When Arnoldi iteration is performed on the Krylov subspace generated using the matrix $\mathbf{A}=\displaystyle \begin{bmatrix}  2& 1& 1& 0\\ 1 &3 &1& 0\\ 0& 1& 3& 1\\ 0& 1& 1& 2 \end{bmatrix}$, the results can depend strongly on the initial vector $\mathbf{u}$. 

    **(a)** Compare the results for `Q` and `H` using the following initial guesses. 
    
    *(i)* `u=[1,0,0,0]`, $\qquad$ *(ii)* `u=[1,1,1,1]`, $\qquad$ *(iii)* `u=rand(4)`.
   
    **(b)** Can you explain what is different about case (ii) in part (a)?  You may have to show the output of {numref}`Function {number}<function-arnoldi>` to see what is happening. 

    (problem-modifiedgs)=
6. ✍ As mentioned in the text, {numref}`Function {number}<function-arnoldi>` does not compute $H_{ij}$ as defined by {eq}`arnoldiip`, but rather 
  
    $$
    S_{ij} = \mathbf{q}_i^* ( \mathbf{A}\mathbf{q}_j - S_{1j}\,\mathbf{q}_1 - \cdots -
    S_{i-1,j}\,\mathbf{q}_{i-1} )
    $$

    for $i=1,\ldots,j$. Show that $S_{ij}=H_{ij}$. (Hence the function is mathematically equivalent to our Arnoldi formulas.)

    ```{index} eigenvalue
    ```
    (problem-arnoldieig)=
7. One way to approximate the eigenvalue problem $\mathbf{A}\mathbf{x}=\lambda\mathbf{x}$ over $\mathcal{K}_m$ is to restrict $\mathbf{x}$ to $\mathcal{K}_m$, so that $\mathbf{x}=\mathbf{Q}_m\mathbf{z}$. Then multiply the approximate equality $\mathbf{A}\mQ_m\mathbf{z} \approx \lambda \mathbf{Q}_m\mathbf{z}$ on the left by $\mathbf{Q}_m^*$. 
 
    **(a)** ✍ Show starting from {eq}`arnoldimat` that

    $$
    \mathbf{Q}_m^* \mathbf{A} \mathbf{Q}_m =  \tilde{\mathbf{H}}_m,
    $$

    where $\tilde{\mathbf{H}}_m$ is the upper Hessenberg matrix resulting from deleting the last row of $\mathbf{H}_m$. What is the size of this matrix? 

    **(b)** ✍ Show the reasoning above leads to the approximate eigenvalue problem $\tilde{\mathbf{H}}_m\mathbf{z} \approx \lambda\mathbf{z}$. (That is, the exact eigenvalues of $\tilde{\mathbf{H}}_m$ ought to approximate some of the eigenvalues of $\mathbf{A}$.)  

    **(c)** ⌨ Using {numref}`Function {number}<function-arnoldi>` and the matrix of {prf:ref}`demos-subspace-unstable`, compute eigenvalues of $\tilde{\mathbf{H}}_m$ for $m=1,\ldots,40$, keeping track of the error between the largest of those values and the largest eigenvalue of $\mathbf{A}$. Make a semi-log graph of the error as a function of $m$.  
    

