# MINRES and conjugate gradients

```{index} matrix; symmetric
```

We have seen before that certain matrix properties enhance solutions to linear algebra problems. One of the most important of these is when $\mathbf{A}^*=\mathbf{A}$; i.e., $\mathbf{A}$ is hermitian. The Arnoldi iteration has a particularly useful specialization to this case. Starting from {eq}`arnoldimat`, we left-multiply by $\mathbf{Q}_m^*$ to get

:::{math}
\mathbf{Q}_m^* \mathbf{A} \mathbf{Q}_m = \mathbf{Q}_m^* \mathbf{Q}_{m+1} \mathbf{H}_m = \tilde{\mathbf{H}}_m,
:::

where $\tilde{\mathbf{H}}_m$ is rows 1 through $m$ of $\mathbf{H}_m$. If $\mathbf{A}$ is hermitian, then so is the left side of this equation, hence $\tilde{\mathbf{H}}_m$ is hermitian too. But it is also upper Hessenberg, so it must be lower Hessenberg as well. The conclusion is that $\tilde{\mathbf{H}}_m$ is tridiagonal.

```{index} Arnoldi iteration
```
Equation {eq}`arnoldivec` of the Arnoldi iteration now simplifies to a much shorter expression:

:::{math}
:label: lanczos
\mathbf{A} \mathbf{q}_m = H_{m-1,m} \,\mathbf{q}_{m-1} + H_{mm} \,\mathbf{q}_m + H_{m+1,m}\,\mathbf{q}_{m+1}.
:::

```{index} Lanczos iteration
```

As before in deriving the Arnoldi iteration, when given the first $m$ vectors we can solve for the entries in column $m$ of $\mathbf{H}$ and then for $\mathbf{q}_{m+1}$. The resulting process is known as the {term}`Lanczos iteration`. In a simple implementation, it needs just a single minor change from Arnoldi, but numerical stability requires some extra effort. We do not present the details. The most important practical difference is that while Arnoldi needs $O(m)$ steps to get $\mathbf{q}_{m+1}$ from the previous vectors, Lanczos needs only $O(1)$ steps, so restarting isn't required for symmetric matrices.

## MINRES

```{index} GMRES; relation to MINRES
```
```{index} MINRES
```

When $\mathbf{A}$ is hermitian and the Arnoldi iteration is reduced to Lanczos, the analog of GMRES is known as MINRES. Like GMRES, MINRES minimizes the residual $\|\mathbf{b}-\mathbf{A}\mathbf{x}\|$ over increasingly large Krylov spaces.

MINRES is also more theoretically tractable than GMRES. Recall that the eigenvalues of a hermitian matrix are real. Of the eigenvalues that are positive, let $\kappa_+$ be the ratio of the one farthest from the origin (largest) to the one closest to the origin (smallest). Similarly, let $\kappa_-$ be the ratio of the negative eigenvalue farthest from the origin to the negative eigenvalue closest to the origin. Then there is a rigorous upper bound on the residual:

:::{math}
:label: minres-conv
\frac{\|\mathbf{r}_m\|_2}{\|\mathbf{b}\|_2} \le  \left( \frac{\sqrt{\kappa_+\kappa_-} - 1}{\sqrt{\kappa_+\kappa_-} + 1} \right)^{\lfloor m/2\rfloor},
:::

```{index} convergence rate; linear
```
where $\lfloor m/2\rfloor$ means to round $m/2$ down to the nearest integer. This bound (though not necessarily MINRES itself) obeys a linear convergence rate. As the product $\kappa_+\kappa_-$ grows, the rate of this convergence approaches one; i.e., is slower. Hence the convergence of MINRES may depend strongly on the eigenvalues of the matrix, with eigenvalues close to the origin (relative to the max eigenvalues) predicted to force a slower convergence.

## Conjugate gradients

```{index} conjugate gradients
```
```{index} matrix; positive definite
```

Given another property in addition to symmetry, we arrive at perhaps the most famous Krylov subspace method for $\mathbf{A}\mathbf{x}=\mathbf{b}$, called {term}`conjugate gradients`. Suppose now that $\mathbf{A}$ is hermitian and positive definite (HPD). Then $\mathbf{A}$ has a Cholesky factorization, which in the complex case is $\mathbf{A}=\mathbf{R}^*\mathbf{R}$. Therefore, for any vector $\mathbf{u}$,

$$
\mathbf{u}^*\mathbf{A}\mathbf{u} = (\mathbf{R}\mathbf{u})^*(\mathbf{R}\mathbf{u})=\|\mathbf{R} \mathbf{u}\|^2,
$$

which is nonnegative and zero only when $\mathbf{u}=\boldsymbol{0}$, provided $\mathbf{A}$ (and therefore $\mathbf{R}$) is nonsingular. Hence we can define a special vector norm relative to $\mathbf{A}$:

:::{math}
:label: Anorm
\| \mathbf{u} \|_{\mathbf{A}} = \left( \mathbf{u}^*\mathbf{A}\mathbf{u} \right)^{1/2}.
:::

The conjugate gradients algorithm minimizes the error, as measured in the $\mathbf{A}$-norm, over the sequence of Krylov subspaces. That is, $\mathbf{x}_m$ makes $\|\mathbf{x}_m-\mathbf{x}\|_{\mathbf{A}}$ as small as possible over $\mathcal{K}_m$. We do not show any details or code for the resulting algorithm.

## Convergence

```{index} matrix!condition number
```

The convergence of CG and MINRES is dependent on the eigenvalues of $\mathbf{A}$. In the HPD case the eigenvalues are real and positive, and they equal the singular values. Hence the condition number $\kappa$ is equal to the ratio of the largest eigenvalue to the smallest one. The following theorem suggests that MINRES and CG are not so different in convergence.

(theorem-cg-converge)= 
::::{proof:theorem} MINRES and CG convergence
Let $\mathbf{A}$ be real and SPD with 2-norm condition number $\kappa$. For MINRES define $R(m)=\|\mathbf{r}_m\|_2/\|\mathbf{b}\|_2$, and for CG define $R(m)=\|\mathbf{x}_m-\mathbf{x}\|_{\mathbf{A}}/\|\mathbf{x}\|_{\mathbf{A}}$,
where $\mathbf{r}_m$ and $\mathbf{x}_m$ are the residual and solution approximation associated with the space $\mathcal{K}_m$. Then

:::{math}
:label: cgconv
R(m) \le  2 \left( \frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1} \right)^m.
:::
::::

As in the indefinite case with MINRES, a larger eigenvalue ratio (condition number) is associated with slower convergence in the positive definite case for both MINRES and CG.} That is, large condition number is a double penalty, increasing both the time it takes to obtain a solution and the effects of numerical errors. Specifically, to make the bound of {eq}`cgconv` less than a number $\epsilon$ requires

\begin{gather*}
  2 \left( \frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1} \right)^m \approx
  \epsilon \\
  m \log \left( \frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1} \right)
  \approx \log\Bigl( \frac{\epsilon}{2} \Bigr).
\end{gather*}

We estimate

\begin{align*}
   \frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1}
 &=  (1 - \kappa^{-1/2}\,) (1 + \kappa^{-1/2}\,)^{-1}\\
 &= (1 - \kappa^{-1/2}\,)  (1 - \kappa^{-1/2} + \kappa^{-1} + \cdots)\\
 &= 1 - 2\kappa^{-1/2} + O(\kappa^{-1}), \quad \text{ as $\kappa
   \rightarrow \infty$.}
\end{align*}

With the Taylor expansion $\log(1+x) = x - (x^2/2) + \cdots$, we
finally conclude

\begin{gather*}
  2 m \kappa^{-1/2} \approx \log\Bigl( \frac{\epsilon}{2} \Bigr),
  \text{ or }
  m = O(\sqrt{\kappa}),
\end{gather*}

as an estimate of the number of iterations needed to achieve a fixed
accuracy. This estimate fails for very large $\kappa$, however.

::::{prf:example} Julia demo
:class: demo
{doc}`demos/minrescg-conv`
::::

We can explain some of the behavior seen in [our demo](`demos/minrescg-conv`). The first matrix has a condition number of $10^2$, whereas the second has $\kappa=10^4$. The linear convergence bounds of the two cases after 100 iterations have values $(9/11)^{100}\approx 2\times 10^{-9}$ and $(99/101)^{100}\approx 0.14$, respectively, which agrees fairly well with the observed reductions in the residual norms. The major practical difference between MINRES and CG lies in the interpretation of minimization of the residual versus minimization of the error in the $\mathbf{A}$-norm.

<!-- 
\begin{exercises}
	\input{krylov/exercises/Symmetry}
\end{exercises} -->


