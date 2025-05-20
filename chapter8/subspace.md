---
numbering:
  enumerator: 8.4.%s
---
(section-krylov-subspace)=
# Krylov subspaces

The power and inverse iterations have a flaw that seems obvious once it is pointed out. Given a seed vector $\mathbf{u}$, they produce a sequence of vectors $\mathbf{u}_1,\mathbf{u}_2,\ldots$ that are scalar multiples of $\mathbf{u},\mathbf{A}\mathbf{u},\mathbf{A}^{2}\mathbf{u},\ldots$, but only the most recent vector is used to produce an eigenvector estimate. 

It stands to reason that we could do no worse, and perhaps much better, if we searched among all linear combinations of the vectors seen in the past. In other words, we seek a solution in the range (column space) of the matrix

:::{math}
:label: krylovmatrix
\mathbf{K}_m =
\begin{bmatrix}
  \mathbf{u} & \mathbf{A}\mathbf{u} & \mathbf{A}^{2} \mathbf{u} & \cdots & \mathbf{A}^{m-1} \mathbf{u}
\end{bmatrix}.
:::

```{index} ! Krylov matrix, ! Krylov subspace
```

::::{prf:definition} Krylov matrix and subspace
Given $n\times n$ matrix $\mathbf{A}$ and $n$-vector $\mathbf{u}$, the $m$th **Krylov matrix** is the $n\times m$ matrix {eq}`krylovmatrix`. The range (i.e., column space) of this matrix is the $m$th **Krylov subspace** $\mathcal{K}_m$.
::::

In general, we expect that the dimension of the Krylov[^kreeluv] subspace $\mathcal{K}_m$, which is the rank of $\mathbf{K}_m$, equals $m$, though it may be smaller.

[^kreeluv]: The proper pronunciation of "Krylov" is something like "kree-luv," but American English speakers often say "kreye-lahv." 

## Properties

As we have seen with the power iteration, part of the appeal of the Krylov matrix is that it can be generated in a way that fully exploits the sparsity of $\mathbf{A}$, simply through repeated matrix-vector multiplication. Furthermore, we have some important mathematical properties.

(theorem-subspace-krylovmult)=
::::{prf:theorem}
Suppose $\mathbf{A}$ is $n\times n$, $0<m<n$, and a vector $\mathbf{u}$ is used to generate Krylov subspaces. If $\mathbf{x}\in\mathcal{K}_m$, then the following hold:

1. $\mathbf{x} = \mathbf{K}_m \mathbf{z}$ for some $\mathbf{z}\in\mathbb{C}^m$.
2. $\mathbf{x} \in \mathcal{K}_{m+1}$.
3. $\mathbf{A}\mathbf{x} \in \mathcal{K}_{m+1}$.
::::

::::{prf:proof}
:enumerated: false


If $\mathbf{x}\in\mathcal{K}_m$, then for some coefficients $c_1,\ldots,c_m$,

:::{math}
\mathbf{x} = c_1 \mathbf{u} + c_2 \mathbf{A} \mathbf{u} + \cdots + c_m \mathbf{A}^{m-1} \mathbf{u}.
:::

Thus let $\mathbf{z}= \begin{bmatrix} c_1 & \cdots & c_m \end{bmatrix}^T$. Also, $\mathbf{x}\in\mathcal{K}_{m+1}$, as we can add zero times $\mathbf{A}^{m}\mathbf{u}$ to the sum. Finally,
  
:::{math}
\mathbf{A}\mathbf{x} = c_1 \mathbf{A} \mathbf{u} + c_2 \mathbf{A}^{2} \mathbf{u} + \cdots + c_m \mathbf{A}^{m} \mathbf{u} \in \mathcal{K}_{m+1}.
:::
::::

## Dimension reduction

```{index} dimension reduction
```

The problems $\mathbf{A}\mathbf{x}=\mathbf{b}$ and $\mathbf{A}\mathbf{x}=\lambda\mathbf{x}$ are posed in a very high-dimensional space $\mathbb{R}^n$ or $\mathbb{C}^n$. One way to approximate them is to replace the full $n$-dimensional space with a much lower-dimensional $\mathcal{K}_m$ for $m\ll n$. This is the essence of the Krylov subspace approach.

For instance, we can interpret $\mathbf{A}\mathbf{x}_m\approx \mathbf{b}$ in the sense of linear least-squares—that is, using @theorem-subspace-krylovmult to let $\mathbf{x}=\mathbf{K}_m\mathbf{z}$,

:::{math}
:label: gmresdef
\min_{\mathbf{x}\in\mathcal{K}_m} \|  \mathbf{A}\mathbf{x}-\mathbf{b} \|
= \min_{\mathbf{z}\in\mathbb{C}^m} \| \mathbf{A}(\mathbf{K}_m\mathbf{z})-\mathbf{b} \|
= \min_{\mathbf{z}\in\mathbb{C}^m} \| (\mathbf{A}\mathbf{K}_m)\mathbf{z}-\mathbf{b} \|.
:::

The natural seed vector for $\mathcal{K}_m$ in this case is the vector $\mathbf{b}$. In the next example we try to implement {eq}`gmresdef`. We do take one precaution: because the vectors $\mathbf{A}^{k}\mathbf{b}$ may become very large or small in norm, we normalize after each multiplication by $\mathbf{A}$, just as we did in the power iteration.

(demo-subspace-unstable)=
::::{prf:example} Conditioning of the Krylov matrix
`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-subspace-unstable-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-subspace-unstable-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-subspace-unstable-python
:::
````
`````
::::


## The Arnoldi iteration

The breakdown of convergence in {numref}`Demo %s <demo-subspace-unstable>` is due to a critical numerical defect in our approach: the columns of the Krylov matrix {eq}`krylovmatrix` increasingly become parallel to the dominant eigenvector, as {eq}`poweriterconverge` predicts, and therefore to one another. As we saw in {numref}`section-leastsq-qr`, near-parallel vectors create the potential for numerical cancellation. This manifests as a large condition number for $\mathbf{K}_m$ as $m$ grows, eventually creating excessive error when solving the least-squares system.

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

Then the vectors $\mathbf{q}_1,\ldots,\mathbf{q}_m$ are the orthonormal basis we seek for $\mathcal{K}_m$. By @theorem-subspace-krylovmult, we know that $\mathbf{A}\mathbf{q}_m \in \mathcal{K}_{m+1}$, and therefore

:::{math}
:label: arnoldivec
\mathbf{A} \mathbf{q}_m = H_{1m} \, \mathbf{q}_1 + H_{2m} \, \mathbf{q}_2 + \cdots + H_{m+1,m}\, \mathbf{q}_{m+1}
:::

for some choice of the $H_{ij}$. Note that by using orthonormality, we have

:::{math}
:label: arnoldiip
\mathbf{q}_i^* (\mathbf{A}\mathbf{q}_m) = H_{im}, \qquad i=1,\ldots, m.
:::

Since we started by assuming that we know $\mathbf{q}_1,\ldots,\mathbf{q}_m$, the only unknowns in {eq}`arnoldivec` are $H_{m+1,m}$ and $\mathbf{q}_{m+1}$. But they appear only as a product, and we know that $\mathbf{q}_{m+1}$ is a *unit* vector, so they are uniquely defined (up to sign) by the other terms in the equation.

We can now proceed iteratively. 

```{index} ! Arnoldi iteration
```
(algorithm-subspace-arnoldi)=
::::{prf:algorithm} Arnoldi iteration
Given matrix $\mathbf{A}$ and vector $\mathbf{u}$:

1. Let $\mathbf{q}_1= \mathbf{u} \,/\, \| \mathbf{u}\|$.
2. For $m=1,2,\ldots$
    
    a. Use {eq}`arnoldiip` to find $H_{im}$ for $i=1,\ldots,m$.
    
    b. Let

    :::{math}
    :label: arnoldigs
    \mathbf{v} = (\mathbf{A} \mathbf{q}_m) - H_{1m} \, \mathbf{q}_1 - H_{2m}\, \mathbf{q}_2 - \cdots - H_{mm}\, \mathbf{q}_m.
    :::

    c. Let $H_{m+1,m}=\|\mathbf{v}\|$.
    
    d. Let $\mathbf{q}_{m+1}=\mathbf{v}\,/\,H_{m+1,m}$.

Return the $n\times (m+1)$ matrix $\mathbf{Q}_{m+1}$ whose columns are $\mathbf{q}_1,\dots, \mathbf{q}_{m+1}$ and the $(m+1)\times m$ matrix $\mathbf{H}_m$.
::::

The vectors $\mathbf{q}_1,\dots, \mathbf{q}_m$ are an orthonormal basis for the space $\mathcal{K}_m$, which makes them ideal for numerical computations in that space. The matrix $\mathbf{H}_m$ plays an important role, too, and has a particular "triangular plus" structure,

:::{math}
:label: arnoldiH
\mathbf{H}_m = \begin{bmatrix}
    H_{11} & H_{12} & \cdots & H_{1m} \\
    H_{21} & H_{22} & \cdots & H_{2m} \\
    & H_{32} & \ddots & \vdots \\
    & & \ddots & H_{mm} \\
    & & & H_{m+1,m}
\end{bmatrix}.
:::

```{index} ! upper Hessenberg matrix
```

::::{prf:definition} Upper Hessenberg matrix
A matrix $\mathbf{H}$ is an {term}`upper Hessenberg matrix` if $H_{ij}=0$ whenever $i>j+1$.
::::

The identity @arnoldivec used over all the iterations can be collected into a single matrix equation, 

:::{math}
:label: arnoldimat
\mathbf{A}\mathbf{Q}_m = \mathbf{Q}_{m+1} \mathbf{H}_m,
:::

which is a fundamental identity of Krylov subspace methods.

## Implementation

An implementation of the Arnoldi iteration is given in {numref}`Function {number} <function-arnoldi>`. A careful inspection shows that inner nested loop does not exactly implement {eq}`arnoldiip` and {eq}`arnoldigs`. The reason is numerical stability. Though the described and implemented versions are mathematically equivalent in exact arithmetic (see @problem-subspace-modifiedgs), the approach in {numref}`Function {number} <function-arnoldi>` is more stable.

(function-arnoldi)=
``````{prf:algorithm} arnoldi
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #function-arnoldi-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-arnoldi-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #function-arnoldi-python
:::
````
`````
``````

(demo-subspace-arnoldi)=
::::{prf:example} Arnoldi iteration
`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-subspace-arnoldi-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-subspace-arnoldi-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-subspace-arnoldi-python
:::
````
`````
::::

In the next section, we revisit the idea of approximately solving $\mathbf{A}\mathbf{x}=\mathbf{b}$ over a Krylov subspace as suggested in @demo-subspace-arnoldi. A related idea explored in @problem-subspace-arnoldieig is used to approximate the eigenvalue problem for $\mathbf{A}$, which is the approach that underlies `eigs` for sparse matrices.

## Exercises

``````{exercise}
:label: problem-subspace-permute
✍ Let $\mathbf{A}=\displaystyle \begin{bmatrix}
0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \\ 1 & 0 & 0 & 0
\end{bmatrix}.$

**(a)** Find the Krylov matrix $\mathbf{K}_3$ for the seed vector $\mathbf{u}=\mathbf{e}_1$. 

**(b)** Find $\mathbf{K}_3$ for the seed vector $\mathbf{u}=\begin{bmatrix}1; \: 1;\: 1; \: 1\end{bmatrix}.$  
``````

``````{exercise}
:label: problem-subspace-conditioning
⌨ For each matrix, make a table of the 2-norm condition numbers $\kappa(\mathbf{K}_m)$ for $m=1,\ldots,10$. Use a vector of all ones as the Krylov seed.

**(a)** Matrix from {numref}`Demo %s <demo-subspace-unstable>`

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
\end{bmatrix} \: (200\times 200) $
``````

% must stay as #3

``````{exercise}
:label: problem-subspace-matrixpolykrylov
✍ Show that if $\mathbf{x}\in\mathcal{K}_m$, then $\mathbf{x}=p(\mathbf{A})\mathbf{u}$ for a polynomial $p$ of degree at most $m-1$. (See {eq}`matrixpoly` for applying a polynomial to a matrix.)
``````

``````{exercise}
:label: problem-subspace-flops
✍ Compute the asymptotic flop requirements for {numref}`Function {number} <function-arnoldi>`. Assume that due to sparsity, a matrix-vector multiplication $\mathbf{A}\mathbf{u}$ requires only $c n$ flops for a constant $c$, rather than the usual $O(n^2)$. 
``````

``````{exercise}
:label: problem-subspace-initial
⌨ When Arnoldi iteration is performed on the Krylov subspace generated using the matrix $\mathbf{A}=\displaystyle \begin{bmatrix}  2& 1& 1& 0\\ 1 &3 &1& 0\\ 0& 1& 3& 1\\ 0& 1& 1& 2 \end{bmatrix}$, the results can depend strongly on the initial vector $\mathbf{u}$. 

**(a)** Apply {numref}`Function {number} <function-arnoldi>` for 3 iterations and output `Q` (which should be square) and `H` when using the following seed vectors. 

*(i)* $\bigl[1,\,0,\,0,\,0\bigr]$ $\qquad$ *(ii)* $\bigl[1,\,1,\,1,\,1\bigr]$ $\qquad$ *(iii)* `rand(4)`

**(b)** Can you explain why case (ii) in part (a) cannot finish successfully? (Hint: What line(s) of the function can possibly return NaN when applied to finite values?) 
``````

% must stay as #6

``````{exercise}
:label: problem-subspace-modifiedgs
✍ As mentioned in the text, {numref}`Function {number} <function-arnoldi>` does not compute $H_{ij}$ as defined by {eq}`arnoldiip`, but rather 

$$
S_{ij} = \mathbf{q}_i^* ( \mathbf{A}\mathbf{q}_j - S_{1j}\,\mathbf{q}_1 - \cdots -
S_{i-1,j}\,\mathbf{q}_{i-1} )
$$

for $i=1,\ldots,j$. Show that $S_{ij}=H_{ij}$. (Hence the function is mathematically equivalent to our Arnoldi formulas.)
``````

```{index} eigenvalue
```

% must stay as #7

``````{exercise}
:label: problem-subspace-arnoldieig
One way to approximate the eigenvalue problem $\mathbf{A}\mathbf{x}=\lambda\mathbf{x}$ over $\mathcal{K}_m$ is to restrict $\mathbf{x}$ to the low-dimensional spaces $\mathcal{K}_m$. 

**(a)** ✍ Show starting from {eq}`arnoldimat` that

$$
\mathbf{Q}_m^* \mathbf{A} \mathbf{Q}_m =  \tilde{\mathbf{H}}_m,
$$

where $\tilde{\mathbf{H}}_m$ is the upper Hessenberg matrix resulting from deleting the last row of $\mathbf{H}_m$. What is the size of this matrix? 

**(b)** ✍ Show the reasoning above leads to the approximate eigenvalue problem $\tilde{\mathbf{H}}_m\mathbf{z} \approx \lambda\mathbf{z}$. (Hint: Start with $\mathbf{A}\mathbf{x} \approx \lambda\mathbf{x}$, and let $\mathbf{x}=\mathbf{Q}_m\mathbf{z}$ before applying part (a).)  

**(c)** ⌨ Apply {numref}`Function {number} <function-arnoldi>` to the matrix of {numref}`Demo %s <demo-subspace-unstable>` using a random seed vector. Compute eigenvalues of $\tilde{\mathbf{H}}_m$ for $m=1,\ldots,40$, keeping track in each case of the error between the largest of those values (in magnitude) and the largest eigenvalue of $\mathbf{A}$. Make a log-linear graph of the error as a function of $m$.  
``````
