---
numbering:
  enumerator: 8.8.%s
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
---
```{code-cell}
:tags: [remove-cell]
clear all
format short
set(0, 'defaultaxesfontsize', 12)
set(0, 'defaultlinelinewidth', 1.5)
set(0, 'defaultFunctionLinelinewidth', 1.5)
set(0, 'defaultscattermarkerfacecolor', 'flat')
gcf;
set(gcf, 'Position', [0 0 600 350])
addpath FNC-matlab
```

(section-krylov-precond)=

# Preconditioning

An important aspect of MINRES and CG (and, by extension, GMRES) is that the convergence of a Krylov method can be expected to deteriorate as the condition number of the matrix increases. Even moderately large condition numbers can make the convergence impractically slow. Therefore, it's common for these methods to be used with a technique to reduce the relevant condition number.

```{index} GMRES; preconditioning in, ! preconditioning
```

::::{prf:definition} Preconditioner
:label: definition-preconditioner
Given a linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$, a {term}`preconditioner` is a matrix $\mathbf{M}$ or equivalent linear transformation that modifies the system to be

```{math}
:label: precond
(\mathbf{M}^{-1} \mathbf{A}) \mathbf{x} = \mathbf{M}^{-1}\mathbf{b}.
```

::::

More specifically, {eq}`precond` is known as *left preconditioning*, which is the simplest and most common type.

As usual, we do not want to actually compute $\mathbf{M}^{-1}$ for a given $\mathbf{M}$. Instead, we have a linear system with the matrix $\mathbf{M}^{-1}\mathbf{A}$. In a Krylov method, the operation "let $\mathbf{v}=\mathbf{A}\mathbf{u}$" becomes a two-step process:

1. Set $\mathbf{y}=\mathbf{A}\mathbf{u}$.
2. Solve $\mathbf{M}\mathbf{v}=\mathbf{y}$ for $\mathbf{v}$.

```{index} sparse matrix, LU factorization
```

As an implementation detail, it is common to provide the Krylov solver with code that does step 2; if the matrix $\mathbf{M}$ is given, the default is to use sparse factorization.

There are competing objectives in the choice of $\mathbf{M}$. On one hand, we want $\mathbf{M}^{-1}\mathbf{A}\approx \mathbf{I}$ in some sense because that makes {eq}`precond` easy to solve by Krylov iteration. Hence, $\mathbf{M}\approx \mathbf{A}$. On the other hand, we desire that solving the system $\mathbf{M}\mathbf{v}=\mathbf{y}$ be relatively fast.

:::{prf:observation}
Good preconditioning is a matter of finding an easily inverted (i.e., quickly solvable) approximation of the original matrix.
:::

## Diagonal preconditioning

One of the simplest choices for the preconditioner $\mathbf{M}$ is a diagonal matrix. This definitely meets the requirement of being fast to invert: the solution of $\mathbf{M}\mathbf{v}=\mathbf{y}$ is just $v_i=y_i/M_{ii}$. The only question is whether it can be chosen in such a way that $\mathbf{M}^{-1}\mathbf{A}$ is much more amenable to Krylov iterations than $\mathbf{A}$ is. This may be the case when the rows of $\mathbf{A}$ differ greatly in scale, or when $\mathbf{A}$ is diagonally dominant (see {eq}`diag-dominant`).

::::{prf:example} Diagonal preconditioning
:label: demo-precond-diagonal

Here is an SPD matrix that arises from solving partial differential equations.

```{code-cell}
A = gallery("wathen", 60, 60);
n = size(A, 1);
clf,  spy(A)
```

There is an easy way to use the diagonal elements of $\mathbf{A}$, or any other vector, as a diagonal preconditioner.

```{code-cell}
M = spdiags(diag(A), 0, n, n);
```

We now compare MINRES with and without the preconditioner.

```{code-cell}
:tags: hide-input

b = ones(n, 1);
[x, ~, ~, ~, resid_plain] = minres(A, b, 1e-10, 400);
clf,  semilogy(resid_plain)
xlabel('iteration number'), ylabel('residual norm')
title('Unpreconditioned MINRES')

[x, ~, ~, ~, resid_prec] = minres(A, b, 1e-10, 400, M);
hold on,  semilogy(resid_prec)
title('Precondtioned MINRES')
legend('no prec.', 'with prec.');
```

The diagonal preconditioner cut down substantially on the number of iterations. The effect on the total time is less dramatic, but this is not a large version of the problem.


::::

## Incomplete factorization

```{index} ! LU factorization; incomplete
```

Another general-purpose technique is the **incomplete LU factorization**. Since true factorization of a sparse matrix usually leads to an undesirable amount of fill-in, incomplete LU sacrifices exact factors by dropping elements smaller than an adjustable threshold.

::::{prf:example} Incomplete LU preconditioning
:label: demo-precond-gmres

Here is a random nonsymmetric matrix.

```{code-cell}
n = 8000;
A = speye(n) + sprand(n, n, 0.00035);
```

Without a preconditioner, restarted GMRES makes slow progress.

```{code-cell}
b = rand(n, 1);
[x, ~, ~, ~, resid_plain] = gmres(A, b, 50, 1e-10, 3);  % restart at 50
format short e
resid_plain(1:30:end)
```

```{index} ! MATLAB; ilu
```

This version of incomplete LU factorization simply prohibits fill-in for the factors, freezing the sparsity pattern of the approximate factors to match the original matrix.

```{code-cell}
[L, U] = ilu(A);
clf
subplot(121), spy(L)
title('L')
subplot(122), spy(U)
title('U')
disp(sprintf("There are %d nonzeros in A", nnz(A)))
```

It does _not_ produce a true factorization of $\mathbf{A}$.

```{code-cell}
norm( full(A - L * U) )  
```

The actual preconditioning matrix is $\mathbf{M}=\mathbf{L}\mathbf{U}$. However, the `gmres` function allows setting the preconditioner by giving the factors independently.

```{code-cell}
[x, ~, ~, ~, resid_prec] = gmres(A, b, [], 1e-10, 300, L, U);
```

The preconditioner makes a significant difference in the number of iterations needed.

```{code-cell}
clf, semilogy(resid_plain)
hold on, semilogy(resid_prec)
xlabel('iteration number'), ylabel('residual norm')
title('Precondtioned GMRES ')
legend('no preconditioner', 'with preconditioner');
```

::::

In practice, good preconditioning is often as important, if not more important, than the specific choice of Krylov method. Effective preconditioning may require deep understanding of the underlying application, however, which limits our ability to go into further details. For instance, the linear system may be some approximation of a continuous mathematical model, and then $\mathbf{M}$ can be derived by using a cruder form of the approximation. Krylov methods offer a natural way to exploit these and other approximate inverses.

## Exercises

``````{exercise}
:label: problem-precond-spd
✍ Suppose $\mathbf{M}=\mathbf{R}^T\mathbf{R}$. Show that the eigenvalues of $\mathbf{R}^{-T}\mathbf{A}\mathbf{R}^{-1}$ are the same as the eigenvalues of $\mathbf{M}^{-1}\mathbf{A}$. (This observation underlies preconditioning variants for SPD matrices.)
``````

``````{exercise}
:label: problem-precond-ilu
⌨ The object returned by `ilu` stores the factors in a way that optimizes sparse triangular substitution. You can recover the factors themselves via

```julia
iLU = ilu(A,τ=0.1)   # for example
L, U = I+iLU.L, iLU.U'
```

In this problem, use `A = 1.5I + sprand(800,800,0.005)`.

**(a)** Using $\tau=0.3$ for the factorization, plot the eigenvalues of $\mathbf{A}$ and of $\mathbf{M}^{-1}\mathbf{A}$ in the complex plane on side-by-side subplots. Do they support the notion that $\mathbf{M}^{-1}\mathbf{A}$ is "more like" an identity matrix than $\mathbf{A}$ is? (Hint: the matrices are small enough to convert to standard dense form for the use of `eigvals`.)

**(b)** Repeat part (a) for $\tau=0.03$. Is $\mathbf{M}$ more accurate than in part (a), or less?
``````

``````{exercise}
:label: problem-precond-surround
⌨ (Continuation of @problem-gmres-surround.) Let $\mathbf{B}$ be `diagm(1:100)`,  let $\mathbf{I}$ be `I(100)`, and let $\mathbf{Z}$ be a $100\times 100$ matrix of zeros. Define 

$$
\mathbf{A} = \begin{bmatrix}
\mathbf{B} & \mathbf{I} \\ \mathbf{Z} & -\mathbf{B}
\end{bmatrix}
$$ 

and let $\mathbf{b}$ be a 200-vector of ones. The matrix $\mathbf{A}$ is difficult for GMRES. 

**(a)** Design a diagonal preconditioner $\mathbf{M}$, with all diagonal elements equal to $1$ or $-1$, such that $\mathbf{M}^{-1}\mathbf{A}$ has all positive eigenvalues. Apply `gmres` without restarts using this preconditioner and a tolerance of $10^{-10}$ for 100 iterations. Plot the convergence curve. 

**(b)** Now design another diagonal preconditioner such that all the eigenvalues of $\mathbf{M}^{-1}\mathbf{A}$ are $1$, and apply preconditioned `gmres` again. How many iterations are apparently needed for convergence? 
``````

``````{exercise}
:label: problem-precond-bai
⌨ Let `A = matrixdepot("Bai/rdb2048")`, and let `b` be a vector of 2048 ones. In the steps below, use GMRES for up to 300 iterations without restarts and with a stopping tolerance of $10^{-4}$.

**(a)** Time the GMRES solution without preconditioning. Verify that convergence was achieved. 

**(b)** Show that diagonal preconditioning is not helpful for this problem.

**(c)** To two digits, find a value of $\tau$ in iLU such that the preconditioned method transitions from effective and faster than part (a) to ineffective. 

``````
