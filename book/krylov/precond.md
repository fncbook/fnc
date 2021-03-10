# Preconditioning

```{index} GMRES; preconditioning
```
```{index} preconditioning
```

An important aspect of MINRES and CG (and, by extension, GMRES) is that the convergence of a Krylov method can be expected to deteriorate as the condition number of the matrix increases. Even moderately large condition numbers can make the convergence impractically slow. Therefore it's common for these methods to be used with a technique known as {term}`preconditioning` to reduce the relevant condition number.

A problem $\mathbf{A}\mathbf{x}=\mathbf{b}$ with a difficult $\mathbf{A}$ can be made more tractable in the mathematically equivalent form

:::{math}
:label: precond
(\mathbf{M}^{-1} \mathbf{A}) \mathbf{x} = \mathbf{M}^{-1}\mathbf{b}
:::

for a matrix $\mathbf{M}$ of our choosing. One goal in this choice is to make $\mathbf{M}^{-1}\mathbf{A}\approx \mathbf{I}$, which makes {eq}`precond` easy to solve by Krylov iteration. In a vague sense, this means $\mathbf{M}\approx \mathbf{A}$. 

On the other hand, there is an important constraint on $\mathbf{M}$. As usual, we do not wish to actually compute $\mathbf{M}^{-1}$. Instead, we have a linear system with the matrix $\mathbf{M}^{-1}\mathbf{A}$, and we take a two-step process to compute any $\mathbf{y}=\mathbf{M}^{-1}\mathbf{A}\mathbf{v}$ within the Krylov iteration:

1. Set $\mathbf{u}=\mathbf{A}\mathbf{v}$.
2. Solve $\mathbf{M}\mathbf{y}=\mathbf{u}$ for $\mathbf{y}$.

Hence we desire that solving the system $\mathbf{M}\mathbf{y}=\mathbf{u}$ be relatively fast. In short, preconditioning is a matter of looking for an inexpensive — that is, easily inverted — approximation of the original matrix.

```{index} matrix; factorization; LU
```
```{index} sparse matrix
```

Methods for deriving a good preconditioner are numerous and often problem-dependent. Certain generic algebraic tricks are available. One of these is an **incomplete LU factorization**. Since true factorization of a sparse matrix usually leads to an undesirable amount of fill-in, incomplete LU prohibits or limits the fill-in in exchange for not getting an exact factorization.

::::{prf:example} Julia demo
:class: demo
:label: demos-precond-gmres
{doc}`demos/precond-gmres`
::::

In many applications the problem $\mathbf{A}\mathbf{x}=\mathbf{b}$ has a known structure that may be exploited. It may be some approximation of a continuous mathematical model, and then $\mathbf{M}$ can be derived by using a cruder form of the approximation. Another important idea is to distinguish near-field and far-field influences in a physically motivated problem and make a fast approximation of the far field. 

Preconditioning is quite important to the practical use of Krylov iterations. However, they usually require deep understanding of the underlying problem that produces the particular linear system to be solved, and we cannot go further into the details here.

## Exercises

(problem-precond-spd)=
1. not available
   
     <!-- ✍ Show that the matrix $\mathbf{R}^{-T}\mathbf{A}\mR^{-1}$ in {eq}`precond-sym` is SPD, given that $\mathbf{A}$ is SPD and $\mathbf{R}$ is nonsingular. -->

2. ✍ Suppose $\mathbf{M}=\mathbf{R}^T\mathbf{R}$. Show that the eigenvalues of $\mathbf{R}^{-T}\mathbf{A}\mathbf{R}^{-1}$ are the same as the eigenvalues of $\mathbf{M}^{-1}\mathbf{A}$.

3. ⌨ Using the definitions in {prf:ref}`demos-precond-gmres`, with $\mathbf{M}=\mathbf{L}\mathbf{U}$, plot the eigenvalues of $\mathbf{A}$ and of $\mathbf{M}^{-1}\mathbf{A}$ in the complex plane. Do they support the notion that $\mathbf{M}^{-1}\mathbf{A}$ is "more like" an identity matrix than $\mathbf{A}$ is? (You have to convert a sparse matrix to dense form using `Matrix(A)` in order to apply `eigvals` to it.)

    ::::{only} solutions
    A = 0.6*speye(1000) + sprand(1000,1000,0.005,1/10000);
    plot(eig(full(A)),'x')

    %%
    [L,U] = ilu(A);
    B = U\(L\A);
    hold on, plot(eig(full(B)),'o')
    % eigenvalues of B are centered at and clustered near 1
    ::::

4. ⌨ (Continuation of [an earlier exercise](problem-gmres-surround).) Let $\mathbf{B}$ be `diagm(1:100)`,  let $\mathbf{I}$ be `I(100)`, and let $\mathbf{Z}$ be a $100\times 100$ matrix of zeros. Define 
  
    $$
    \mathbf{A} = \begin{bmatrix}
      \mathbf{B} & \mathbf{I} \\ \mathbf{Z} & -\mathbf{B}
    \end{bmatrix}
    $$ 
  
    and let $\mathbf{b}$ be a $200\times 1$ vector of ones. The matrix $\mathbf{A}$ is difficult for GMRES. 
  
    **(a)** Design a diagonal preconditioner $\mathbf{M}$ such that $\mathbf{M}^{-1}\mathbf{A}$ has all positive eigenvalues. Apply `gmres` without restarts using this preconditioner and a tolerance of $10^{-10}$ for 100 iterations. Plot the convergence curve. 
  
    **(b)** Now design another diagonal preconditioner such that all the eigenvalues of $\mathbf{M}^{-1}\mathbf{A}$ are 1 and apply preconditioned `gmres` again. How many iterations are apparently needed for full convergence? 

    ::::{only} solutions
    B = diag(1:100);
    A = [B eye(100); zeros(100) -B];
    b = ones(200,1);
    M = diag([ones(100,1);-ones(100,1)]);
    [~,~,~,~,rv] = gmres(A,b,100,1e-10,1,M);
    semilogy(rv)
    %% (b)
    M = diag([(1:100),-(1:100)]);
    [~,~,~,~,rv] = gmres(A,b,100,1e-10,1,M);
    semilogy(rv)
    % depending on how you count, it's 2 or 3 iterations. 
    ::::

