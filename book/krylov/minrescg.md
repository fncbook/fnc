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

::::{prf:theorem} MINRES and CG convergence
:label: theorem-cg-converge
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
:label: demos-minrescg-conv
{doc}`demos/minrescg-conv`
::::

We can explain some of the behavior seen in [our demo](`demos/minrescg-conv`). The first matrix has a condition number of $10^2$, whereas the second has $\kappa=10^4$. The linear convergence bounds of the two cases after 100 iterations have values $(9/11)^{100}\approx 2\times 10^{-9}$ and $(99/101)^{100}\approx 0.14$, respectively, which agrees fairly well with the observed reductions in the residual norms. The major practical difference between MINRES and CG lies in the interpretation of minimization of the residual versus minimization of the error in the $\mathbf{A}$-norm.

## Exercises

1. ✍ For each part, the eigenvalues of $\mathbf{A}$ are given. Suppose MINRES is applied to solve $\mathbf{A}\mathbf{x}=\mathbf{b}$. Find a numerical value for the upper bound in {eq}`minres-conv` or {eq}`cgconv`, whichever is most appropriate. Then determine which of the cases gives the slowest and which gives the fastest convergence.

    **(a)** $-100,-99,\ldots,-1,1,2,\ldots,100$
    
    **(b)** $-100,1,2,\ldots,100$
    
    **(c)** $1,2,\ldots,100$

    ::::{only} solutions
    In (a) you have $\kappa_-=\kappa_+=100$, and the bound is $(99/101)^(m/2)$. In (b) you have $\kappa_-=1$ and $\kappa_+=100$, so the bound is $(9/10)^(m/2)$. In (c) it's SPD with $\kappa=100$ and the bound is $(9/10)^m=(.81)^(m/2)$. So (a) is slowest and (c) is fastest.
    ::::

2. ⌨ Let $\mathbf{b}$ be a random unit vector of length 200. Define the matrix
   
    ```julia
    u = LinRange(-200,-1,100); 
    v = LinRange(10,100,100);
    A = diagm([u;v]);
    ```
    
    **(a)**  Apply 120 iterations of `minres` to solve $\mathbf{A}\mathbf{x}=\mathbf{b}$. Compute the relative norm of the answer. Plot the norm of the residual as a function of $m$.

    **(b)** Add to your graph the line representing the upper bound {eq}`minres-conv`. (Ignore the rounding in the exponent.) This line should stay strictly on or above the convergence curve.

    ::::{only} solutions
    u = linspace(-200,-1); v = linspace(10,100);
    A = diag([u v]);
    b = rand(200,1);  b = b/norm(b);
    [x,~,~,~,rr] = minres(A,b,[],120);
    clf, semilogy(0:120,rr)
    err = norm(A\b-x)/norm(A\b)
    hold on
    km = 200; kp = 10;
    K = (sqrt(km*kp)-1)/(sqrt(km*kp)+1);
    semilogy(K.^((0:120)/2),'–')
    ::::

3. ⌨ Let $\mathbf{b}$ be a random unit vector of length 500. Define the matrix
   
    ```julia
    A = spdiagm(LinRange(4,10000,500));
    ```

    **(a)**  Apply 80 iterations of `minres` to solve $\mathbf{A}\mathbf{x}=\mathbf{b}$. Compute the relative norm of the answer. Plot the norm of the residual as a function of $m$.
    
    **(b)** Add to your graph the line representing the upper bound {eq}`cgconv`. This line should stay strictly on or above the convergence curve.
  
    **(c)** Add a convergence curve for 80 iterations of `cg`.

    ::::{only} solutions
    A = diag(linspace(4,10000,500));
    b = rand(500,1);  b = b/norm(b);
    [x,~,~,~,rr] = minres(A,b,[],80);
    clf, semilogy(0:80,rr)
    norm(A\b-x)/norm(A\b)
    hold on
    k = 2500;
    K = (sqrt(k)-1)/(sqrt(k)+1);
    semilogy(2*K.^((0:80)),'–')
    [x,~,~,~,rr] = pcg(A,b,[],80);
    semilogy(0:80,rr)
    ::::

4. ✍ Given real $n\times n$ symmetric $\mathbf{A}$ and vector $\mathbf{b}=\mathbf{A}\mathbf{x}$, we can define the scalar-valued function
  
    :::{math}
    \varphi(\mathbf{u}) = \mathbf{u}^T \mathbf{A} \mathbf{u} - 2 \mathbf{u}^T \mathbf{b}, \qquad \mathbf{u}\in\mathbb{R}^n.
    :::

    **(a)** Expand and simplify the expression $\varphi(\mathbf{x}+\mathbf{v})-\varphi(\mathbf{x})$, keeping in mind that $\mathbf{A}\mathbf{x}=\mathbf{b}$.
  
    **(b)** Using the result of (a), prove that if $\mathbf{A}$ is an SPD matrix, $\varphi$ has a global minimum at $\mathbf{x}$.
  
    **(c)** Show that for any vector $\mathbf{u}$, $\|\mathbf{u}-\mathbf{x}\|_{\mathbf{A}}^2-\varphi(\mathbf{u})$ is constant.
  
    **(d)** Using the result of (c), prove that CG minimizes $\varphi(\mathbf{u})$ over Krylov subspaces.

5. ⌨ The following linear system arises from the Helmholtz equation for wave propagation:
    
    ```julia
    A = FNC.poisson(n) - k^2*I;
    b = -ones(n^2);
    ```

    **(a)** Repeat {prf:ref}`demos-minrescg-conv` using this linear system with $n=50$ and $k=1.3\pi$.
  
    **(b)** Repeat {prf:ref}`demos-minrescg-conv` using this linear system with $n=50$ and $k=1.5\pi$. Use "eig" to explain why CG fails in this case.

    ::::{only} solutions
    n = 50;
    k = 1.3*pi;
    A = n^2*gallery('poisson',n) - k^2*speye(n^2);
    b = -ones(n^2,1);
    [xMR,~,~,~,residMR] = minres(A,b,1e-12,100);
    [xCG,~,~,~,residCG] = pcg(A,b,1e-12,100);
    clf
    semilogy(0:100,residMR/norm(b),'.-')
    hold on
    semilogy(0:length(residCG)-1,residCG/norm(b),'.-')
    ::::


