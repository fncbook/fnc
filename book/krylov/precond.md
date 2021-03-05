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

::::{proof:example} Julia demo
:class: demo
{doc}`demos/precond-gmres`
::::

In many applications the problem $\mathbf{A}\mathbf{x}=\mathbf{b}$ has a known structure that may be exploited. It may be some approximation of a continuous mathematical model, and then $\mathbf{M}$ can be derived by using a cruder form of the approximation. Another important idea is to distinguish near-field and far-field influences in a physically motivated problem and make a fast approximation of the far field. 

Preconditioning is quite important to the practical use of Krylov iterations. However, they usually require deep understanding of the underlying problem that produces the particular linear system to be solved, and we cannot go further into the details here.

<!-- 

\begin{exercises}
  \input{krylov/exercises/Preconditioning}
\end{exercises} -->

