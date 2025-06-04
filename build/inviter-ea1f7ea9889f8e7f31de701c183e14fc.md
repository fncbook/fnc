---
numbering:
  enumerator: 8.3.%s
---
(section-krylov-inviter)=

# Inverse iteration

Power iteration finds only the dominant eigenvalue. We next show that it can be adapted to find any eigenvalue, provided you start with a reasonably good estimate of it. Some simple linear algebra is all that is needed.

::::{prf:theorem}
Let $\mathbf{A}$ be an $n\times n$ matrix with eigenvalues $\lambda_1,\ldots,\lambda_n$ (possibly with repeats), and let $s$ be a complex scalar. Then:

1. The eigenvalues of the matrix $\mathbf{A}-s\mathbf{I}$ are $\lambda_1-s,\ldots,\lambda_n-s$.
2. If $s$ is not an eigenvalue of $\mathbf{A}$, the eigenvalues of the matrix $(\mathbf{A}-s\mathbf{I})^{-1}$ are $(\lambda_1-s)^{-1},\ldots,(\lambda_n-s)^{-1}$.
3. The eigenvectors associated with the eigenvalues in the first two parts are the same as those of $\mathbf{A}$.
::::

::::{prf:proof}
:enumerated: false

The equation $\mathbf{A}\mathbf{v}=\lambda \mathbf{v}$ implies that $(\mathbf{A}-s\mathbf{I})\mathbf{v} = \mathbf{A}\mathbf{v} - s\mathbf{I}\mathbf{v} = \lambda\mathbf{v} - s\mathbf{v} = (\lambda-s)\mathbf{v}$. That proves the first part of the theorem. For the second part, we note that by assumption, $(\mathbf{A}-s\mathbf{I})$ is nonsingular, so $(\mathbf{A}-s\mathbf{I})\mathbf{v} = (\lambda-s) \mathbf{v}$ implies that $\mathbf{v} = (\lambda-s) (\mathbf{A}-s\mathbf{I}) \mathbf{v}$, or $ (\lambda-s)^{-1} \mathbf{v} =(\mathbf{A}-s\mathbf{I})^{-1} \mathbf{v}$. The discussion above also proves the third part of the theorem.
::::

Consider first part 2 of the theorem with $s=0$, and suppose that $\mathbf{A}$ has a *smallest* eigenvalue,

```{math}
|\lambda_n| \ge |\lambda_{n-1}| \ge \cdots > |\lambda_1|.
```

Then clearly

```{math}
|\lambda_1^{-1}| > |\lambda_{2}^{-1}| \ge \cdots \ge |\lambda_n^{-1}|,
```

```{index} eigenvalue; dominant
```

and $\mathbf{A}^{-1}$ has a {term}`dominant eigenvalue`. Hence, power iteration on $\mathbf{A}^{-1}$ can be used to find the eigenvalue of $\mathbf{A}$ closest to zero. For nonzero values of $s$, then we suppose there is an ordering

```{math}
:label: shiftorder
|\lambda_n-s| \ge \cdots \ge |\lambda_2-s|  > |\lambda_1-s|.
```

Then it follows that

```{math}
|\lambda_1-s|^{-1} > |\lambda_{2}-s|^{-1} \ge \cdots \ge |\lambda_n-s|^{-1},
```

and power iteration on the matrix $(\mathbf{A}-s\mathbf{I})^{-1}$ converges to $(\lambda_1-s)^{-1}$, which is easily solved for $\lambda_1$ itself.

## Algorithm

A literal application of @definition-poweriteration would include the step

```{math}
:label: shiftinvstepbad
\mathbf{y}_k = (\mathbf{A}-s\mathbf{I})^{-1} \mathbf{x}_k.
```

As always, however, we do not want to explicitly find the inverse of a matrix. Instead, we should implement this step as the solution of a linear system.

```{index} ! inverse iteration
```

```{index} see: shifted inverse iteration; inverse iteration
```

::::{prf:definition} Inverse iteration
:label: definition-inviter
Given matrix $\mathbf{A}$ and shift $s$:

1. Choose $\mathbf{x}_1$ such that $\twonorm{\mathbf{x}_1} = 1$.
2. For $k=1,2,\ldots$,

    a. Solve for $\mathbf{y}_k$ in
    :::{math}
    :label: shiftinvstep
    (\mathbf{A}-s\mathbf{I}) \mathbf{y}_k =\mathbf{x}_k.
    :::

    b. Set $\beta_k = s + \dfrac{1}{\mathbf{x}_k^* \mathbf{y}_k}$. 

    c. Set $\mathbf{x}_{k+1} = {\mathbf{y}_k} / {\twonorm{\mathbf{y}_k}}.$

Return $\beta_1,\beta_2,\ldots$ as eigenvalue estimates, and $\mathbf{x}_1,\mathbf{x}_2,\ldots$ as associated eigenvector estimates.
::::

:::{note}
In @definition-poweriteration, we used $\mathbf{x}_k^* \mathbf{y}_k$ as an estimate of the dominant eigenvalue of $\mathbf{A}$. Here, that ratio is an estimate of $(\lambda_1-s)^{-1}$, and solving for $\lambda_1$ gives the $\beta_k$ in @definition-inviter.
:::

Each pass of inverse iteration requires the solution of a linear system of equations with the matrix $\mathbf{B}=\mathbf{A}-s\mathbf{I}$. This solution might use methods we consider later in this chapter. Here, we use (sparse) PLU factorization and hope for the best. Since the matrix $\mathbf{B}$ is constant, the factorization needs to be done only once for all iterations. The details are in {numref}`Function {number} <function-inviter>`.

``````{prf:algorithm} inviter
:label: function-inviter

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #function-inviter-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-inviter-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #function-inviter-python
:::
````
`````
``````

## Convergence rate

The convergence is linear, at a rate found by reinterpreting @theorem-poweriterconv with $(\mathbf{A}-s\mathbf{I})^{-1}$ in place of $\mathbf{A}.$ With the eigenvalues ordered as in {eq}`shiftorder`, in the general case we have

```{math}
:label: inviterconv
\frac{\abs{\beta_{k+1} - \lambda_1}}{\abs{\beta_{k} - \lambda_1}} \rightarrow
\abs{\frac{  \lambda_1 - s } {\lambda_2 - s}}\quad \text{ as } \quad k\rightarrow \infty,
```

and in the hermitian case, we have

```{math}
:label: inviterconvherm
\frac{\abs{\beta_{k+1} - \lambda_1}}{\abs{\beta_{k} - \lambda_1}} \rightarrow
\abs{\frac{  \lambda_1 - s } {\lambda_2 - s}}^2  \quad \text{ as } \quad k\rightarrow \infty.
```

Thus, the convergence is best when the shift $s$ is close to the target eigenvalue—specifically, when it is much closer to that eigenvalue than to any other.

::::{prf:example} Convergence of inverse iteration
:label: demo-inviter-conv

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-inviter-conv-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-inviter-conv-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-inviter-conv-python
:::
```` 
`````

::::

## Rayleigh quotient iteration

There is a clear opportunity for positive feedback in @definition-inviter. The convergence rate of inverse iteration improves as the shift gets closer to the true eigenvalue—and the algorithm computes improving eigenvalue estimates! Updating the shift to $s=\beta_k$ after each iteration greatly accelerates the convergence. You are asked to implement this algorithm in @problem-inviter-dynamicshift.

```{index} convergence rate; quadratic
```

If the eigenvalues are ordered by distance to $s$, then (asymptotically) one step of inverse iteration reduces the error by the factor $|\lambda_1-s|/|\lambda_2-s|$. As $s \to\lambda_1$, the change in the denominator is negligible. So, if at one point the error $\abs{\lambda_1-s}$ is about $\epsilon$, then the error in the next estimate is reduced by a factor $O(\epsilon)$, making it $O(\epsilon^2)$. That is, each step now squares the error, which is {term}`quadratic convergence`.

::::{prf:example} Dynamic shift strategy
:label: demo-inviter-accel

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-inviter-accel-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-inviter-accel-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-inviter-accel-python
:::
```` 
`````

::::

There is a price to pay for this improvement. The matrix of the linear system to be solved, $(\mathbf{A}-s\mathbf{I}),$ now changes with each iteration. That means that we can no longer do just one LU factorization for the entire iteration. The speedup in convergence usually makes this tradeoff worthwhile, however.

In practice power and inverse iteration are not as effective as the algorithms used by `eigs` and based on the mathematics described in the rest of this chapter. However, inverse iteration can be useful for turning an eigenvalue estimate into an eigenvector estimate.

## Exercises

``````{exercise}
:label: problem-inviter-convergence
⌨  Use {numref}`Function {number} <function-inviter>` to perform 10 iterations for the given matrix and shift. Compare the results quantitatively to the convergence given by {eq}`inviterconv`.

**(a)**  $\mathbf{A} = \begin{bmatrix}
1.1 & 1 \\
0 & 2.1
\end{bmatrix}, \; s = 1 \qquad $
**(b)** $\mathbf{A} = \begin{bmatrix}
1.1 & 1 \\
0 & 2.1
\end{bmatrix}, \; s = 2\qquad $

**(c)** $\mathbf{A} = \begin{bmatrix}
1.1 & 1 \\
0 & 2.1
\end{bmatrix}, \; s = 1.6\qquad $
**(d)** $\mathbf{A} = \begin{bmatrix}
2 & 1 \\
1 & 0
\end{bmatrix}, \; s = -0.33 \qquad$

**(e)** $\mathbf{A} = \begin{bmatrix}
6 & 5 & 4 \\
5 & 4 & 3 \\
4 & 3 & 2
\end{bmatrix}, \;  s = 0.1 $
``````

``````{exercise}
:label: problem-inviter-shift
✍ Let $\mathbf{A} = \displaystyle \begin{bmatrix} 1.1 & 1 \\ 0 & 2.1 \end{bmatrix}.$ Given the starting vector $\mathbf{x}_1=[1,1]$, find the vector $\mathbf{x}_2$ for the following shifts.

**(a)** $s=1\quad$ **(b)** $s=2\quad$ **(c)** $s=1.6$
``````

``````{exercise}
:label: problem-inviter-imaginary
✍ Why is it a bad idea to use unshifted inverse iteration with the matrix $\displaystyle \begin{bmatrix} 0 & 1 \\ -1 & 0 \end{bmatrix}$? Does the shift $s=-1$ improve matters?
``````

``````{exercise}
:label: problem-inviter-illconditioned
✍ When the shift $s$ is very close to an eigenvalue of $\mathbf{A}$, the matrix $\mathbf{A}-s\mathbf{I}$ is close to a singular matrix. But then {eq}`shiftinvstep` is a linear system with a badly conditioned matrix, which should create a lot of error in the numerical solution for $\mathbf{y}_k$. However, it happens that the error is mostly in the direction of the eigenvector we are looking for, as the following toy example illustrates.

Prove that $\displaystyle \begin{bmatrix} 1 & 1 \\ 0 & 0 \end{bmatrix}$ has an eigenvalue at zero with associated eigenvector $\mathbf{v}=[-1,1]^T$. Suppose this matrix is perturbed slightly to $\displaystyle \mathbf{A} = \begin{bmatrix} 1 & 1 \\ 0 & \epsilon \end{bmatrix}$, and that $\mathbf{x}_k=[1,1]$ in {eq}`shiftinvstep`. Show that once $\mathbf{y}_k$ is normalized by its infinity norm, the result is within $\epsilon$ of a multiple of $\mathbf{v}$.
``````

% must stay as #5

``````{exercise}
:label: problem-inviter-lumpmembrane
⌨ (Continuation of @problem-power-lumpmembraneeig.) This exercise concerns the $n^2\times n^2$ sparse matrix defined by `FNC.poisson(n)` for integer $n$. It represents a lumped model of a vibrating square membrane held fixed around the edges.

**(a)** The eigenvalues of $\mathbf{A}$ closest to zero are approximately squares of the frequencies of vibration for the membrane. Using `eigs`, find the eigenvalue $\lambda_m$ closest to zero for $n=10,15,20,25$.

**(b)** For each $n$ in part (a), apply 50 steps of {numref}`Function {number} <function-inviter>` with zero shift. On one graph, plot the four convergence curves $|\beta_k-\lambda_m|$ using a semi-log scale.

**(c)** Let `v` be the eigenvector (second output) found by {numref}`Function {number} <function-inviter>` for $n=25$. Make a surface plot of the vibration mode by reshaping `v` into an $n\times n$ matrix.
``````

% must remain as number 6

``````{exercise}
:label: problem-inviter-dynamicshift
⌨ This problem explores the use of Rayleigh quotient iteration.

**(a)** Modify {numref}`Function {number} <function-inviter>` to change the value of the shift $s$ to be the most recently computed value in the vector $\beta$. Note that the matrix `B` must also change with each iteration, so the LU factorization cannot be done just once.

**(b)** Define a $100\times 100$ matrix with values $k^2$ for $k=1,\ldots,100$ on the main diagonal and random values uniformly distributed between 0 and 1 on the first superdiagonal. (Since this matrix is triangular, the diagonal values are its eigenvalues.) Using an initial shift of $s=920$, apply Rayleigh quotient iteration. Determine which eigenvalue was found and make a table of the `log10` of the errors in the iteration as a function of iteration number. (These should approximately double, until machine precision is reached, due to quadratic convergence.)

**(c)** Repeat part (b) using a different initial shift of your choice.
``````
