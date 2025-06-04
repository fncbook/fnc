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

Analysis of matrix powers is most straightforward in the diagonalizable case. Let $\mathbf{A}$ be any diagonalizable $n\times n$ matrix having eigenvalues $\lambda_1,\ldots,\lambda_n$ and corresponding linearly independent eigenvectors $\mathbf{v}_1,\ldots,\mathbf{v}_n$. For our later convenience (and without losing any generality), we assume that the eigenvectors are normalized so that

```{math}
\twonorm{\mathbf{v}_j} = 1 \quad \text{for } j=1,\ldots,n.
```

We also make an important assumption about the eigenvalue magnitudes that will hold for most but not all matrices.

::::{prf:definition} Dominant eigenvalue
:label: definition-dominanteigenvalue
If the eigenvalues of a matrix are such that

```{math}
:label: evorder
|\lambda_1| > |\lambda_2| \ge |\lambda_3| \ge \cdots \ge |\lambda_n|,
```

then we say that $\lambda_1$ is the {term}`dominant eigenvalue` of the matrix.
::::

In @demo-power-one, for instance, $\lambda_1=1$ is the dominant eigenvalue.

Now let $\mathbf{x}$ be an $n$-vector, let $k$ be a positive integer, and refer to {eq}`evdpower`:

```{math}
\mathbf{A}^k \mathbf{x} = \mathbf{V}\mathbf{D}^k\mathbf{V}^{-1}\mathbf{x}.
```

Let $\mathbf{z}=\mathbf{V}^{-1}\mathbf{x}$, and recall that $\mathbf{D}$ is a diagonal matrix of eigenvalues. Then

````{math}
:label: powerAkx0
\begin{split}
  \mathbf{A}^k\mathbf{x} &= \mathbf{V}\mathbf{D}^k \mathbf{z} = \mathbf{V}\begin{bmatrix} \lambda_1^kz_1 \\[0.5ex] \lambda_2^kz_2 \\ \vdots \\ \lambda_n^kz_n \end{bmatrix} \\
 &= \lambda_1^k \left[ z_1 \mathbf{v}_{1} +
  z_2 \left(\frac{\lambda_2}{\lambda_1}\right) ^k 
  \mathbf{v}_{2} + \cdots + z_n \left(\frac{\lambda_n}{\lambda_1}\right)^k
  \mathbf{v}_{n} \right].
\end{split}
````

Since $\lambda_1$ is dominant, we conclude that if $z_1\neq 0$,

```{math}
:label: poweriterconverge
\begin{split}
\twonorm{ \frac{ \mathbf{A}^k\mathbf{x}}{z_1 \lambda_1^k} - \mathbf{v}_1 } 
& \le \underbrace{\left|\frac{z_2}{z_1}\right|}_{c_2} \cdot \underbrace{ \left| \frac{\lambda_2}{\lambda_1} \right|^k}_{r_2^k}  \twonorm{\mathbf{v}_{2}}
+ \cdots + \underbrace{\left|\frac{z_n}{z_1}\right|}_{c_n} \cdot \underbrace{\left|\frac{\lambda_n}{\lambda_1}\right|^k}_{r_n^k}  \twonorm{\mathbf{v}_{n}}  \\[1mm]
& = \sum_{j=2}^n c_j r_j^k  \cdot 1 \\ 
& \rightarrow 0 \text{ as $k\rightarrow \infty$},
\end{split}
```

since, by @evorder, $r_j = | \lambda_j / \lambda_1 | < 1$ for $j=2,\ldots,n$. If we choose $\mathbf{x}$ randomly, then the probability that $z_1=0$ is zero, so we will not be concerned with that case.

```{important}
If $\mathbf{A}$ is diagonalizable and has a dominant eigenvalue $\lambda_1$, then $\mathbf{A}^k\mathbf{x}$ asymptotically approaches a scalar multiple of the dominant eigenvector.
```

:::{attention}
For algorithmic purposes, it is important to interpret $\mathbf{A}^k\mathbf{x}$ as $\mathbf{A}\bigl( \cdots\bigl( \mathbf{A} (\mathbf{A}\mathbf{x})\bigl) \cdots\bigl)$, i.e., as repeated applications of $\mathbf{A}$ to a vector. Doing so allows us to fully exploit sparsity of $\mathbf{A}$, something which is not preserved by taking a matrix power $\mathbf{A}^k$ explicitly before the multiplication with $\mathbf{x}$ (see @demo-structure-fill).
:::

## Algorithm

An important technicality separates us from an algorithm: unless $|\lambda_1|=1$, the factor $\lambda_1^k$ tends to make $\|\mathbf{A}^k\mathbf{x}\|$ either very large or very small. In practice, we cannot easily normalize by $\lambda_1^k$ as we did in {eq}`poweriterconverge`, since we don't know $\lambda_1$ in advance.

This issue is resolved by alternating matrix–vector multiplications with renormalizations.

```{index} ! power iteration
```

::::{prf:definition} Power iteration
:label:  definition-poweriteration
Given matrix $\mathbf{A}$:

1. Choose $\mathbf{x}_1$ such that $\twonorm{\mathbf{y}_k}=1$.
2. For $k=1,2,\ldots$,

    a. Set $\mathbf{y}_k = \mathbf{A} \mathbf{x}_k$.

    b. Set $\beta_k = \mathbf{x}_k^* \mathbf{y}_k$. (In the real case, this is $\beta_k = \mathbf{x}_k^T \mathbf{y}_k$.)

    c. Set $\mathbf{x}_{k+1} = \dfrac{\mathbf{y}_k}{\twonorm{\mathbf{y}_k}}$.

Return $\beta_1,\beta_2,\ldots$ as dominant eigenvalue estimates, and $\mathbf{x}_1,\mathbf{x}_2,\ldots$ as associated eigenvector estimates.
::::

:::{tip}
The vectors and scalars in @definition-poweriteration are subscripted by iteration number to help make the discussion here more convenient. In practice, the algorithm can be implemented without keeping the history.
:::

Observe that we can write

```{math}
:label: powernorm
\mathbf{x}_{k} = (\alpha_1 \alpha_2 \cdots \alpha_k ) \mathbf{A}^k \mathbf{x}_{1},
```

where $\alpha_j = \twonorm{\mathbf{y}_j}^{-1}$ is the normalization factor at iteration $j$. Thus, the reasoning of @poweriterconverge implies that $\mathbf{x}_k$ converges to a scalar multiple of the dominant eigenvector $\mathbf{v}_1$ of $\mathbf{A}$. Specifically, suppose that there is some number $\gamma$ with $|\gamma|=1$ such that

```{math}
 \mathbf{x}_k - \gamma \mathbf{v}_1 = \epsilon \mathbf{w}, 
 ```

 where $\norm{\mathbf{w}} = 1$ and $\epsilon \ll 1$. Then

```{math}
:label: eq-powereverror
\begin{split}
\beta_k &= \mathbf{x}_k^* \mathbf{y}_k  \\ 
&= \mathbf{x}_k^* \mathbf{A} \mathbf{x}_k \\ 
&= \bigl(\overline{\gamma} \mathbf{v}_1^* + \epsilon \mathbf{w}^*\bigr) \mathbf{A} \bigl(\gamma \mathbf{v}_1 + \epsilon \mathbf{w}\bigr) \\ 
& = |\gamma|^2 \mathbf{v}_1^* (\lambda_1 \mathbf{v}_1) + O(\epsilon) \\
& = \lambda_1 + O(\epsilon).
\end{split}
```

That is, $\beta_k$ from the power iteration estimates $\lambda_1$ about as well as $\mathbf{x}_k$ estimates the dominant eigenvector.

````{note}
If $\mathbf{A}$ is hermitian, then $\beta_k$ equals the {term}`Rayleigh quotient` $R_{\mathbf{A}}(\mathbf{x}_k)$, and, according to @theorem-rayleighquotient, 

```{math}
:label: eq-powersquare
\beta_k = \lambda_1 + O(\epsilon^2),
```

so the eigenvalue estimate is even better.
````

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

```{important}
The only use of $\mathbf{A}$ is to find the matrix-vector product $\mathbf{A}\mathbf{x}$, which makes exploitation of the sparsity of $\mathbf{A}$ trivial.
```

## Convergence rate

While we now feel confident that the sequence $\{\beta_k\}$ converges to the dominant eigenvalue $\lambda_1$, we would like to know how fast this convergence is. Looking back at @powerAkx0, we know that our normalizations make it so that

```{math}
:label: eq-powererr
\mathbf{x}_{k} - \gamma \mathbf{v}_1 =  b_2 \left(\frac{\lambda_2}{\lambda_1}\right) ^k \mathbf{v}_{2} 
+ \cdots + b_n \left(\frac{\lambda_n}{\lambda_1}\right)^k \mathbf{v}_{n}
```

for some constants $b_2,\ldots,b_n$. If we now make a stronger assumption that $\lambda_2$ dominates the rest of the eigenvalues, i.e.,

```{math}
:label: evorder2
|\lambda_1| > |\lambda_2| > |\lambda_3| \ge \cdots \ge |\lambda_n|,
```

then expression on the right-hand side of {eq}`eq-powererr` is dominated by its first term, because

```{math}
\sum_{j=2}^n b_j \left(\frac{\lambda_j}{\lambda_1}\right) ^k \mathbf{v}_{j} 
= \left(\frac{\lambda_2}{\lambda_1}\right)^k \left[ b_2 \mathbf{v}_2 + \underbrace{\sum_{j=3}^n {b_j} \left(\frac{\lambda_j}{\lambda_2}\right) ^k \mathbf{v}_{j} }_{\to 0\, \text{ as } k\to\infty} \right] .
```

Therefore, @eq-powereverror now implies that $|\beta_k - \lambda_1|$ is dominated by $|\lambda_2 / \lambda_1|^k$, which is a case of {term}`linear convergence`.

```{index} convergence rate; linear
```

::::{prf:theorem} Power iteration convergence
:label: theorem-poweriterconv
Given @evorder2, the eigenvalue estimates $\beta_k$ satisfy 

```{math}
:label: poweriterconv
\frac{\abs{\beta_{k+1} - \lambda_1}}{\abs{\beta_{k}-\lambda_1}} \rightarrow \abs{\frac{\lambda_2}{\lambda_1}} \quad \text{as } k\rightarrow \infty
```

in the general case, and

```{math}
:label: poweriterconvherm
\frac{\abs{\beta_{k+1} - \lambda_1}}{\abs{\beta_{k}-\lambda_1}} \rightarrow \abs{\frac{\lambda_2}{\lambda_1}}^2 \quad \text{as } k\rightarrow \infty
```

if $\mathbf{A}$ is hermitian.
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
