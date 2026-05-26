---
numbering:
  enumerator: 8.8.%s
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---
```{code-cell}
:tags: [remove-cell]
from numpy import *
from scipy import linalg
from scipy.linalg import norm
from matplotlib.pyplot import *
from prettytable import PrettyTable
from timeit import default_timer as timer
import sys
sys.path.append('fncbook/')
import fncbook as FNC

# This (optional) block is for improving the display of plots.
# from IPython.display import set_matplotlib_formats
# set_matplotlib_formats("svg","pdf")
# %config InlineBackend.figure_format = 'svg'
rcParams["figure.figsize"] = [7, 4]
rcParams["lines.linewidth"] = 2
rcParams["lines.markersize"] = 4
rcParams['animation.html'] = "jshtml"  # or try "html5"
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
from scipy.sparse import sparray
import rogues
A = rogues.wathen(60, 60)
n = A.shape[0]
print(f"Matrix is {n} x {n} with {A.nnz} nonzeros")
```

There is an easy way to use the diagonal elements of $\mathbf{A}$, or any other vector, as a diagonal preconditioner.

```{code-cell}
import scipy.sparse as sp
prec = sp.diags(1 / A.diagonal(), 0)
```

We now compare CG with and without the preconditioner.

```{code-cell}
:tags: [hide-input]

from scipy.sparse.linalg import cg
b = ones(n)
hist = lambda x: resid.append(norm(b - A @ x))
resid = [norm(b)]
start = timer()
x, _ = cg(A, b, rtol=1e-4, maxiter=200, callback=hist)
print(f"No preconditioner: Finished in {timer() - start:.2f} sec")
resid_plain = resid.copy()
resid = [norm(b)]
start = timer()
x, _ = cg(A, b, rtol=1e-4, maxiter=200, M=prec, callback=hist)
print(f"Diagonal preconditioner: Finished in {timer() - start:.2f} sec")
resid_prec = resid.copy()

semilogy(resid_plain, label="no preconditioner")
semilogy(resid_prec, label="diagonal preconditioner")
xlabel("iteration"), ylabel("residual norm")
legend(),  title("Convergence of CG with and without preconditioning");
```

The diagonal preconditioner cut down substantially on the number of iterations and the execution time.

::::

## Incomplete factorization

```{index} ! LU factorization; incomplete
```

Another general-purpose technique is the **incomplete LU factorization**. Since true factorization of a sparse matrix usually leads to an undesirable amount of fill-in, incomplete LU sacrifices exact factors by dropping elements smaller than an adjustable threshold.

::::{prf:example} Incomplete LU preconditioning
:label: demo-precond-gmres

Here is a random nonsymmetric matrix.

```{code-cell}
import scipy.sparse as sp
n = 8000
A = 2.8 * sp.eye(n) + sp.rand(n, n, 0.002)
```

Without a preconditioner, GMRES can solve a system with this matrix.

```{code-cell}
from scipy.sparse.linalg import gmres

b = random.rand(n)
hist = lambda rvec: resid.append(norm(rvec))
resid = [1.]

start = timer()
x, flag = gmres(A, b, maxiter=300, rtol=1e-10, restart=50, callback=hist)
print(f"time for plain GMRES: {timer() - start:.3f} sec")
resid_plain = resid.copy()
```

```{index} ! Python; spilu
```

The following version of incomplete LU factorization drops all sufficiently small values (i.e., replaces them with zeros). This keeps the number of nonzeros in the factors under control.

```{code-cell}
from scipy.sparse.linalg import spilu
iLU = spilu(A, drop_tol=0.2)
print(f"Factors have {iLU.nnz} nonzeros, while A has {A.nnz}")
```

The result is not a true factorization of the original matrix. However, it's close enough for an approximate inverse in a preconditioner. 

```{code-cell}
from scipy.sparse.linalg import LinearOperator
prec = LinearOperator((n, n), matvec=lambda y: iLU.solve(y))

resid = [1.];  start = timer()
x, flag = gmres(A, b, M=prec, maxiter=300, rtol=1e-10, restart=50, callback=hist)
print(f"time for preconditioned GMRES: {timer() - start:.3f} sec")
resid_prec = resid
```

```{code-cell}
:tags: [hide-input]
semilogy(resid_plain, label="no prec.")
semilogy(resid_prec, label="iLU prec.")
xlabel("iteration number"),  ylabel("residual norm")
legend()
title("GMRES convergence compared");
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
