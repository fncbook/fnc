# Preconditioning

An important aspect of MINRES and CG (and, by extension, GMRES) is that the convergence of a Krylov method can be expected to deteriorate as the condition number of the matrix increases. Even moderately large condition numbers can make the convergence impractically slow. Therefore it's common for these methods to be used with a technique known as 
```{index} GMRES!preconditioning
```
 
```{index} preconditioning
```
 {term}`preconditioning` to reduce the relevant condition number.

A problem $\mathbf{A}\mathbf{x}=\mathbf{b}$ with a difficult $\mathbf{A}$ can be made more tractable in the mathematically equivalent form

:::{math}
  :label: precond
  (\mathbf{M}^{-1} \mathbf{A}) \mathbf{x} = \mathbf{M}^{-1}\mathbf{b}
:::

for a matrix $\mathbf{M}$ of our choosing. One goal in this choice is to make $\mathbf{M}^{-1}\mathbf{A}\approx \mathbf{I}$, which makes {eq}`precond` easy to solve by Krylov iteration. In a loose sense, this means $\mathbf{M}\approx \mathbf{A}$. On the other hand, there is an important constraint on $\mathbf{M}$. As usual, we do not wish to actually compute $\mathbf{M}^{-1}$. Instead, we have a
linear system with the matrix $\mathbf{M}^{-1}\mathbf{A}$, and we take
a two-step process to compute any $\mathbf{y}=\mathbf{M}^{-1}\mathbf{A}\mathbf{v}$ within the Krylov iteration:
\begin{enumerate}
\item Set $\mathbf{u}=\mathbf{A}\mathbf{v}$.
\item Solve $\mathbf{M}\mathbf{y}=\mathbf{u}$ for $\mathbf{y}$.
\end{enumerate}
Hence we desire that solving the system $\mathbf{M}\mathbf{y}=\mathbf{u}$ be relatively fast. In short, \texthighlight{precond}{preconditioning is a matter of looking for an inexpensive—that is, easily inverted—approximation of the original matrix.}

Methods for deriving a good preconditioner are numerous and often problem-dependent. Certain generic algebraic
tricks are available. One of these is an 
```{index} matrix!factorization!LU
```
 
```{index} sparse matrix
```
 **incomplete LU factorization**. Since true factorization of a sparse matrix usually leads to an undesirable amount of fill-in, incomplete LU prohibits or limits the fill-in in exchange for not getting an exact factorization.


::::{proof:example}
  \inputexample{krylov}{precondilu}
::::

In many applications the problem $\mathbf{A}\mathbf{x}=\mathbf{b}$ has a known structure that may be exploited.
It may be some approximation of a continuous mathematical model, and then $\mathbf{M}$ can be derived by using a cruder form of the approximation. Another important idea is to distinguish near-field and far-field influences in a physically motivated problem and make a fast approximation of the far field. These are advanced topics and require appreciation of the underlying problem that produces the linear system to be solved.
<!-- 

\begin{exercises}
  \input{krylov/exercises/Preconditioning}
\end{exercises} -->

