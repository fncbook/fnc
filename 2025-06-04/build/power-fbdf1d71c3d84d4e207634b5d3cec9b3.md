---
numbering:
  enumerator: 8.2.%s
---
(section-krylov-power)=
# Power iteration

```{index} sparse matrix
```

Given that matrix-vector multiplication is fast for sparse matrices, let's see what we might accomplish with only that at our disposal.

::::{prf:example} Power iteration
:label: demo-power-one
`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-power-one-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-power-one-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-power-one-python
:::
````
`````
::::

There was a little cheating in @demo-power-one to make the story come out neatly (specifically, the normalization step after creating a random matrix). But it illustrates an important general fact that we investigate now.

## Dominant eigenvalue

```{index} ! eigenvalue; dominant
```

Analysis of matrix powers is most straightforward in the diagonalizable case. Let $\mathbf{A}$ be any diagonalizable $n\times n$ matrix having eigenvalues $\lambda_1,\ldots,\lambda_n$ and corresponding linearly independent eigenvectors $\mathbf{v}_1,\ldots,\mathbf{v}_n$. We also make an important assumption about the eigenvalue magnitudes.

::::{prf:definition} Dominant eigenvalue
:label: definition-dominant-eigenvalue
If the eigenvalues of a matrix are such that

:::{math}
:label: evorder
|\lambda_1| > |\lambda_2| \ge |\lambda_3| \ge \cdots \ge |\lambda_n|,
:::

then we say that $\lambda_1$ is the {term}`dominant eigenvalue` of the matrix.
::::

In @demo-power-one, for instance, $\lambda_1=1$ is the dominant eigenvalue.

Now let $\mathbf{x}$ be an $n$-vector, let $k$ be a positive integer, and refer to {eq}`evdpower`: 

:::{math}
\mathbf{A}^k \mathbf{x} = \mathbf{V}\mathbf{D}^k\mathbf{V}^{-1}\mathbf{x}.
:::

Let $\mathbf{z}=\mathbf{V}^{-1}\mathbf{x}$, and recall that $\mathbf{D}$ is a diagonal matrix of eigenvalues. Then

::::{math}
:label: powerAkx0
\begin{split}
  \mathbf{A}^k\mathbf{x} &= \mathbf{V}\mathbf{D}^k \mathbf{z} = \mathbf{V}\begin{bmatrix} \lambda_1^kz_1 \\[0.5ex] \lambda_2^kz_2 \\ \vdots \\ \lambda_n^kz_n \end{bmatrix} \\
	&= \lambda_1^k \left[ z_1 \mathbf{v}_{1} +
		z_2 \left(\frac{\lambda_2}{\lambda_1}\right) ^k 
		\mathbf{v}_{2} + \cdots + z_n \left(\frac{\lambda_n}{\lambda_1}\right)^k
		\mathbf{v}_{n} \right].
\end{split}
::::

Since $\lambda_1$ is dominant, we conclude that if $z_1\neq 0$,

:::{math}
:label: poweriterconverge
\left\| \frac{ \mathbf{A}^k\mathbf{x}}{\lambda_1^k}
- z_1\mathbf{v}_1\right\| \le |z_2|\cdot\left|\frac{\lambda_2}{\lambda_1}\right| ^k
\| \mathbf{v}_{2} \| + \cdots +  |z_n|\cdot\left|\frac{\lambda_n}{\lambda_1}\right|^k
\| \mathbf{v}_{n} \| \rightarrow 0 \text{ as $k\rightarrow \infty$}.
:::

That is, $\mathbf{A}^k\mathbf{x}$ eventually is close to close to a scalar multiple of the dominant eigenvector.[^zeromeasure]

[^zeromeasure]: If $\mathbf{x}$ is chosen randomly, the probability that $z_1=0$ is mathematically zero.

:::{attention}
For algorithmic purposes, it is important to interpret $\mathbf{A}^k\mathbf{x}$ as $\mathbf{A}\bigl( \cdots\bigl( \mathbf{A} (\mathbf{A}\mathbf{x})\bigl) \cdots\bigl)$, i.e., as repeated applications of $\mathbf{A}$ to a vector. Doing so allows us to fully exploit sparsity of $\mathbf{A}$, something which is not preserved by taking a matrix power $\mathbf{A}^k$ explicitly before the multiplication with $\mathbf{x}$ (see @demo-structure-fill).
:::
## Power iteration

An important technicality separates us from an algorithm: unless $|\lambda_1|=1$, the factor $\lambda_1^k$ tends to make $\|\mathbf{A}^k\mathbf{x}\|$ either very large or very small. Nor can we easily normalize by $\lambda_1^k$, as in {eq}`poweriterconverge`, unless we know $\lambda_1$ in advance.

To make a practical algorithm, we alternate matrix-vector multiplication with a renormalization of the vector. In the following, we use $x_{k,m}$ and $y_{k,m}$ to mean the $m$th component of vectors $\mathbf{x}_k$ and $\mathbf{y}_k$.

```{index} ! power iteration
```

::::{prf:algorithm} Power iteration
:label: algorithm-power-power
Given matrix $\mathbf{A}$:

1. Choose $\mathbf{x}_1$.
2. For $k=1,2,\ldots$, 
   
    a. Set $\mathbf{y}_k = \mathbf{A} \mathbf{x}_k$.
	
    b. Find $m$ such that $|y_{k,m}|=\|{\mathbf{y}_k} \|_\infty$.

    c. Set $\alpha_k = \dfrac{1}{y_{k,m}}$ and $\,\beta_k = \dfrac{y_{k,m}}{x_{k,m}}$.

    d. Set $\mathbf{x}_{k+1} = \alpha_k \mathbf{y}_k$.

Return $\beta_1,\beta_2,\ldots$ as dominant eigenvalue estimates, and $\mathbf{x}_1,\mathbf{x}_2,\ldots$ as associated eigenvector estimates.
::::

:::{note}
The vectors and scalars in @algorithm-power-power are subscripted by iteration number to help make the discussion here more convenient. In practice, the algorithm can be implemented without keeping the history, and $\alpha_k$ need not be separately computed at all.
:::

By construction, $\| \mathbf{x}_{k}\|_\infty=1$ for all $k > 1$. Also, we can write

:::{math}
:label: powernorm
\mathbf{x}_{k} = (\alpha_1 \alpha_2 \cdots \alpha_k ) \mathbf{A}^k \mathbf{x}_{1}.
:::

Thus {numref}`Algorithm {number} <algorithm-power-power>` modifies {eq}`powerAkx0` and {eq}`poweriterconverge` only slightly.

Finally, if $\mathbf{x}_k$ is nearly a dominant eigenvector of $\mathbf{A}$, then $\mathbf{A}\mathbf{x}_k$ is nearly $\lambda_1\mathbf{x}_k$, and we can take the ratio $\beta_k=y_{k,m}/x_{k,m}$ as an eigenvalue estimate. In fact, revisiting {eq}`powerAkx0`, the extra $\alpha_j$ normalization factors cancel in the ratio, and, after some simplification, we get

:::{math}
:label: poweriterratio
\beta_k = \frac{y_{k,m}}{x_{k,m}} = \lambda_1
\frac{1+r_2^{k+1} b_2 + \cdots +  r_n^{k+1} b_n}{1+r_2^{k} b_2 +  \cdots +  r_n^{k} b_n},
:::

where $r_j=\lambda_j/\lambda_1$ and the $b_j$ are constants. By assumption {eq}`evorder`, each $r_j$ satisfies $|r_j|<1$, so we see that $\beta_k\rightarrow \lambda_1$ as $k\rightarrow\infty$.

{numref}`Function {number} <function-poweriter>` is our implementation of power iteration.

``````{prf:algorithm} poweriter
:label: function-poweriter
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #function-poweriter-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-poweriter-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #function-poweriter-python
:::
````
`````
``````

Observe that the only use of $\mathbf{A}$ is to find the matrix-vector product $\mathbf{A}\mathbf{x}$, which makes exploitation of the sparsity of $\mathbf{A}$ trivial.

## Convergence

Let's examine the terms in the numerator and denominator of {eq}`poweriterratio` more carefully:

:::{math}
:label: powerleftover
\begin{split}
r_2^{k} b_2 +  \cdots +  r_n^{k} b_n &= r_2^k \left[ b_2 + \left( \frac{r_3}{r_2} \right)^kb_3 + \cdots + \left( \frac{r_n}{r_2} \right)^kb_n \right] \\
&= r_2^k \left[ b_2 + \left( \frac{\lambda_3}{\lambda_2} \right)^kb_3 + \cdots + \left( \frac{\lambda_n}{\lambda_2} \right)^kb_n \right].
\end{split}
:::

At this point we'll introduce an additional assumption,

:::{math}
:label: evorder2
|\lambda_2| > |\lambda_3| \ge \cdots \ge |\lambda_n|.
:::

This condition isn't strictly necessary, but it simplifies the following statements considerably because now it's clear that the quantity in {eq}`powerleftover` approaches $b_2 r_2^k$ as $k\rightarrow \infty$.

Next we estimate {eq}`poweriterratio` for large $k$, using a geometric series expansion for the denominator to get

:::{math}
:label: powerest
\begin{split}
\beta_k & \to \lambda_1 \left( 1+b_2 r_2^{k+1} \right) \left( 1 - b_2 r_2^{k} + O(r_2^{2k}) \right), \\
\beta_k - \lambda_1 &\to \lambda_1 b_2 (  r_2 - 1 ) r_2^{k}.
\end{split}
:::

```{index} convergence rate; linear
```
This is {term}`linear convergence` with factor $r_2$:

:::{math}
:label: poweriterconv
\frac{\beta_{k+1} - \lambda_1}{\beta_{k}-\lambda_1} \rightarrow r_2 = \frac{\lambda_2}{\lambda_1} \quad \text{as } k\rightarrow \infty.
:::

::::{prf:observation}
The error in the power iteration eigenvalue estimates $\beta_k$ is reduced asymptotically by a constant factor $\lambda_2/\lambda_1$ at each iteration, where $\lambda_1$ and $\lambda_2$ are the dominant eigenvalues of $\mathbf{A}$.
::::

::::{prf:example} Convergence of power iteration
:label: demo-power-iter
`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-power-iter-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-power-iter-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-power-iter-python
:::
````
`````
::::


The practical utility of {eq}`poweriterconv` is limited: if we knew $\lambda_1$ and $\lambda_2$, we wouldn't be running the power iteration in the first place! Sometimes it's possible to find estimates of or bounds on the ratio. If nothing else, though, it is useful to know that linear convergence is expected at a rate based solely on the dominant eigenvalues. 

## Exercises

``````{exercise}
:label: problem-power-convergence
⌨ Use {numref}`Function {number} <function-poweriter>` to perform 20 power iterations for the following matrices. Quantitatively compare the observed convergence to the prediction in {eq}`poweriterconv`.

**(a)**
$\mathbf{A} = \begin{bmatrix}
1.1 & 1 \\
0.1 & 2.4
\end{bmatrix} \quad$
**(b)** $\mathbf{A} = \begin{bmatrix}
2 & 1 \\
1 & 0
\end{bmatrix} \quad$
**(c)** $ \mathbf{A} = \begin{bmatrix}
6 & 5 & 4 \\
5 & 4 & 3 \\
4 & 3 & 2
\end{bmatrix}$

**(d)** $\mathbf{A} = \begin{bmatrix}
8 & -14 & 0 & -14 \\
-8 & 1 & 1 & 1 \\
-4 & -2 & 0 & 2 \\
8 & -7 & -1 & -7 
\end{bmatrix}$
``````

``````{exercise}
:label: problem-power-stuck
✍ Describe what happens during power iteration using the matrix $\mathbf{A}= \begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix}$ and initial vector $\mathbf{x}=\begin{bmatrix} 0.4\\0.7 \end{bmatrix}$. Does the algorithm converge to an eigenvector? How does this relate to {eq}`powerAkx0`?
``````

``````{exercise}
:label: problem-power-lumpmembraneeig
⌨  In  @problem-linearsystems-lumpstring we considered a mass-lumped model of a hanging string that led to a tridiagonal system of linear equations. Then, in @problem-evd-lumpstring, we found that eigenvectors of the same matrix correspond to vibrational modes of the string. The same setup can be applied to a membrane hanging from a square frame. Lumping the mass onto a Cartesian grid, each interacts with the four neighbors to the north, south, east, and west. If $n$ masses are used in each coordinate direction, we get an $n^2\times n^2$ sparse matrix $\mathbf{A}$ that can be constructed by `FNC.poisson(n)`.

**(a)** Let $n=10$ and make a `spy` plot of $\mathbf{A}$. What is the density of $\mathbf{A}$? Most rows all have the same number of nonzeros; find this number.

**(b)** Find the dominant $\lambda_1$ using `eigs` for $n=10,15,20,25$.

**(c)** For each $n$ in part (b), apply 100 steps of {numref}`Function {number} <function-poweriter>`. On one graph, plot the four convergence curves $|\beta_k-\lambda_1|$ using a semi-log scale. (They will not be smooth curves because this matrix has many repeated eigenvalues that complicate the convergence analysis.) 
``````

``````{exercise}
:label: problem-power-actors
⌨ Copy the instructions from @problem-structure-actorsmat to obtain a large, sparse matrix $\mathbf{A}$. Use {numref}`Function {number} <function-poweriter>` to find the leading eigenvalue of $\mathbf{A}^T\mathbf{A}$ to at least six significant digits.
``````

``````{exercise}
:label: problem-power-rayleigh
⌨ For symmetric matrices, the Rayleigh quotient {eq}`rayleigh` converts an $O(\epsilon)$ eigenvector estimate into an $O(\epsilon^2)$ eigenvalue estimate. Duplicate {numref}`Function {number} <function-poweriter>` and rename it to `powersym`. Modify the new function to use the Rayleigh quotient to produce the entries of `β` or `beta`. Your function should not introduce any additional matrix-vector multiplications. Apply the original {numref}`Function {number} <function-poweriter>` and the new `powersym` on the `MatrixDepot`/`gallery`/`rogues` matrix `fiedler(100)`, plotting the convergence curves on one graph.
``````
