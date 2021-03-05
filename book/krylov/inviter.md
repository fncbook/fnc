# Inverse iteration

Power iteration finds only the dominant eigenvalue. We next show that it can be adapted to find any eigenvalue, provided you start with a reasonably good estimate of it. Some simple linear algebra is all that is needed.


::::{prf:theorem}
Let $\mathbf{A}$ be an $n\times n$ matrix with eigenvalues $\lambda_1,\ldots,\lambda_n$ (possibly with repeats), and let $s$ be a complex scalar. Then:

1. The eigenvalues of the matrix $\mathbf{A}-s\mathbf{I}$ are $\lambda_1-s,\ldots,\lambda_n-s$.
2. If $s$ is not an eigenvalue of $\mathbf{A}$, the eigenvalues of the matrix $(\mathbf{A}-s\mathbf{I})^{-1}$ are $(\lambda_1-s)^{-1},\ldots,(\lambda_n-s)^{-1}$.
3. The eigenvectors associated with the eigenvalues in the first two parts are the same as those of $\mathbf{A}$.
::::

::::{prf:proof}
The equation $\mathbf{A}\mathbf{v}=\lambda \mathbf{v}$ implies that $(\mathbf{A}-s\mathbf{I})\mathbf{v} = \mathbf{A}\mathbf{v} - s\mathbf{I}\mathbf{v} = \lambda\mathbf{v} - s\mathbf{v} = (\lambda-s)\mathbf{v}$. That proves the first part of the theorem. For the second part, we note that by assumption, $(\mathbf{A}-s\mathbf{I})$ is nonsingular, so $(\mathbf{A}-s\mathbf{I})\mathbf{v} = (\lambda-s) \mathbf{v}$ implies that $\mathbf{v} = (\lambda-s) (\mathbf{A}-s\mathbf{I}) \mathbf{v}$, or $ (\lambda-s)^{-1} \mathbf{v} =(\mathbf{A}-s\mathbf{I})^{-1} \mathbf{v}$. The discussion above also proves the third part of the theorem.
::::


Consider first part 2 of the theorem with $s=0$, and suppose that $\mathbf{A}$ has a *smallest* eigenvalue,

:::{math}
|\lambda_n| \ge |\lambda_{n-1}| \ge \cdots > |\lambda_1|.
:::

Then clearly

:::{math}
|\lambda_1^{-1}| > |\lambda_{2}^{-1}| \ge \cdots \ge |\lambda_n^{-1}|,
:::

```{index} eigenvalue; dominant
```
```{index} inverse iteration
```

and $\mathbf{A}^{-1}$ has a dominant eigenvalue. Hence power iteration on $\mathbf{A}^{-1}$ can be used to find the eigenvalue of $\mathbf{A}$ closest to zero. This is called  {term}`inverse iteration`. Comparing to {eq}`evorder` and {eq}`poweriterconv`, it is clear that the linear convergence rate of inverse iteration is the ratio

:::{math}
\frac{\lambda_{2}^{-1}}{\lambda_{1}^{-1}} = \frac{\lambda_{1}}{\lambda_{2}}.
:::

These observations generalize easily to nonzero values of $s$. Specifically, if we suppose that the eigenvalues are ordered by their distance to $s$, i.e.,

:::{math}
:label: shiftorder
|\lambda_n-s| \ge \cdots \ge |\lambda_2-s|  > |\lambda_1-s|,
:::

then it follows that

:::{math}
|\lambda_1-s|^{-1} > |\lambda_{2}-s|^{-1} \ge \cdots \ge |\lambda_n-s|^{-1}.
:::

Hence power iteration on the matrix $(\mathbf{A}-s\mathbf{I})^{-1}$ converges to $(\lambda_1-s)^{-1}$. This is sometimes known as *shifted inverse iteration*, though we use "inverse iteration" to refer to the shifted variety too.

## Algorithm

Shifted inverse iteration introduces a major new algorithmic wrinkle. The key step to consider is the matrix-vector multiplication

:::{math}
:label: shiftinvstepbad
\mathbf{y}_k = (\mathbf{A}-s\mathbf{I})^{-1} \mathbf{x}_k.
:::

As always, we do not want to explicitly find the inverse of a matrix. Instead we should write this step as

:::{math}
:label: shiftinvstep
\text{Solve }  (\mathbf{A}-s\mathbf{I}) \mathbf{y}_k =\mathbf{x}_k \text{ for } \mathbf{y}_k.
:::

Each step of inverse iteration therefore requires the solution of a linear system of equations with the matrix $\mathbf{B}=\mathbf{A}-s\mathbf{I}$. The discussion becomes awkward for us at this moment, because the solution of the large, sparse linear system $\mathbf{B}\mathbf{y}=\mathbf{x}$ is something we consider later on in this chapter. In order to get things working, here we use (sparse) PLU factorization for this system and hope for the best. Note that the matrix of the linear system is constant, so the factorization needs to be done only once for all iterations, with only triangular solves being done repeatedly. The details are in the following implementation.

(function-inviter)=
````{proof:function} inviter
**Shifted inverse iteration for the closest eigenvalue.**

```{code-block} julia
:lineno-start: 1
"""
inviter(A,s,numiter)

Perform `numiter` inverse iterations with the matrix `A` and shift
`s`, starting from a random vector, and return a vector of
eigenvalue estimates and the final eigenvector approximation.
"""
function inviter(A,s,numiter)
    n = size(A,1)
    x = normalize(randn(n),Inf)
    gamma = zeros(numiter)
    fact = lu(A - s*I)
    for k = 1:numiter
      y = fact\x
      normy,m = findmax(abs.(y))
      gamma[k] = x[m]/y[m] + s
      x = y/y[m]
    end

    return gamma,x
end
```
````

One additional detail is worth mentioning. Suppose the power iteration is applied to the matrix $(\mathbf{A}-s\mathbf{I})^{-1}$, producing estimates $\beta_k$. These are actually converging to $(\lambda_1-s)^{-1}$, so to recover $\lambda_1$ we compute $\gamma_k = s+\beta_k^{-1}.$

## Convergence

The convergence rate is found by interpreting {eq}`poweriterconv` from the power iteration in the new context:

:::{math}
:label: inviterconv
\frac{\gamma_{k+1} - \lambda_1}{\gamma_{k} - \lambda_1} \rightarrow
\frac{  \lambda_1 - s } {\lambda_2 - s}\quad \text{ as } \quad k\rightarrow \infty.
:::

We observe that the convergence is best when the shift $s$ is close to the target eigenvalue $\lambda_1$, specifically when it is much closer to that eigenvalue than to any other.

::::{prf:example} Julia demo
:class: demo
{doc}`demos/inviter-conv`
::::

## Dynamic shifting

There is a clear opportunity for positive feedback in {ref}`fun-inviter`. The convergence rate of inverse iteration improves as the shift gets closer to the true eigenvalue — and the output of the algorithm is a sequence of improving eigenvalue estimates! If we update the shift to $s=\gamma_k$ after each iteration, the convergence accelerates. You are asked to implement this algorithm in {ref}`prob-inviter-dynamicshift`.

```{index} convergence rate; quadratic
```
Let's estimate the resulting convergence. If the eigenvalues are ordered by distance to $s$, then the convergence is linear with rate $|\lambda_1-s|/|\lambda_2-s|$. As $s\to\lambda_1$, the change in the denominator is negligible. So if the error $(\lambda_1-s)$ is $\epsilon$, then the error in the next estimate is reduced by a factor $O(\epsilon)$. That is, $\epsilon$ becomes $O(\epsilon^2)$, which is 
 *quadratic* convergence.

::::{sidebar}
:class: demo
{doc}`demos/inviter-accel`
::::

There is a price to pay for this improvement. The matrix of the linear system to be solved, $(\mathbf{A}-s\mathbf{I})\mathbf{y}=\mathbf{x}$, now changes with each iteration. That means that we can no longer do just one LU factorization to do the entire iteration. The speedup in convergence usually makes this tradeoff worthwhile, however.

In practice power and inverse iteration are not as effective as the algorithms used by `eigs` and based on the mathematics described in the rest of this chapter. However, inverse iteration can be useful for turning an eigenvalue estimate into an eigenvector estimate.

<!-- 
\begin{exercises}
    \input{krylov/exercises/InverseIteration}
\end{exercises}
 -->
