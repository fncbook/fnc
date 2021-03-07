# GMRES

The most important use of the Arnoldi iteration is to solve the square linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$. The iteration's generation of an orthogonal basis is the key to fixing the trouble [we encountered before](demos/subspace-unstable). Recall that we replaced the linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$ by the lower-dimensional approximation

:::{math}
\min_{\mathbf{x}\in \mathcal{K}_m} \|  \mathbf{A}\mathbf{x}-\mathbf{b}  \| = \min_{\mathbf{z}\in\mathbb{C}^m} \|   \mathbf{A}\mathbf{K}_m\mathbf{z}-\mathbf{b}  \|,
:::

where $\mathbf{K}_m$ is the Krylov matrix generated using $\mathbf{A}$ and the seed vector $\mathbf{b}$.  This method was unstable due to the poor conditioning of $\mathbf{K}_m$, which is a numerically poor basis of $\mathcal{K}_m$. If we use the columns of $\mathbf{Q}_m$ instead as a basis, then we set $\mathbf{x}=\mathbf{Q}_m\mathbf{z}$ and obtain

:::{math}
:label: gmresproblem
\min_{\mathbf{z}\in\mathbb{C}^m}\, \bigl\| \mathbf{A} \mathbf{Q}_m \mathbf{z}-\mathbf{b}  \bigr\|.
:::

From the fundamental Arnoldi identity {eq}`arnoldimat`, this is equivalent to

:::{math}
:label: gmresproblem1
\min_{\mathbf{z}\in\mathbb{C}^m}\, \bigl\| \mathbf{Q}_{m+1} \mathbf{H}_m\mathbf{z}-\mathbf{b} \bigr\|.
:::

Note that $\mathbf{q}_1$ is a unit multiple of $\mathbf{b}$, so $\mathbf{b} = \|\mathbf{b}\| \mathbf{Q}_{m+1}\mathbf{e}_1$. Thus {eq}`gmresproblem1` becomes

:::{math}
:label: gmresproblem2
\min_{\mathbf{z}\in\mathbb{C}^m}\, \bigl\| \mathbf{Q}_{m+1} (\mathbf{H}_m\mathbf{z}-\|\mathbf{b}\|\mathbf{e}_1) \bigr\|.
:::


The least-squares problems {eq}`gmresproblem`,  {eq}`gmresproblem1`, and {eq}`gmresproblem2` are all $n\times m$. But observe that for any $\mathbf{w}\in\mathbb{C}^{m+1}$,

:::{math}
  \|\mathbf{Q}_{m+1}\mathbf{w}\|^2 = \mathbf{w}^*\mathbf{Q}_{m+1}^*\mathbf{Q}_{m+1}\mathbf{w} = \mathbf{w}^*\mathbf{w} = \|\mathbf{w}\|^2.
:::

The first norm in that equation is on $\mathbb{C}^n$, while the last is on the much smaller space $\mathbb{C}^{m+1}$. Hence the least squares problem {eq}`gmresproblem2` is equivalent to

:::{math}
  :label: gmresproblemsmall
  \min_{\mathbf{z}\in\mathbb{C}^m}\, \bigl\| \mathbf{H}_m\mathbf{z}-\|\mathbf{b}\|\mathbf{e}_1 \bigr\|,
:::

which is of size $(m+1)\times m$. We call the solution of this minimization $\mathbf{z}_m$, and then $\mathbf{x}_m=\mathbf{Q}_m \mathbf{z}_m$ is the $m$th approximation to the solution of $\mathbf{A}\mathbf{x}=\mathbf{b}$.

```{index} GMRES
```
The algorithm resulting from this discussion is known as {term}`GMRES`, for Generalized Minimum RESidual. GMRES uses the output of the Arnoldi iteration to minimize the residual of $\mathbf{A}\mathbf{x}=\mathbf{b}$ over successive Krylov subspaces.


::::{prf:example} Julia demo
:class: demo
:label: demos-gmres-intro
{doc}`demos/gmres-intro`
::::


Compare the graph in [this demo](demos/gmres-intro.ipynb)  to the one in [the earlier unstable version](demos/subspace-unstable.ipynb). Both start with the same linear convergence, but only the version using Arnoldi avoids the instability created by the poor Krylov basis.

A basic implementation of GMRES is given below.

(function-gmres)=
````{proof:function} gmres
**GMRES for a linear system.**

```{code-block} julia
:lineno-start: 1
"""
gmres(A,b,m)

Do `m` iterations of GMRES for the linear system `A`*x=`b`. Return
the final solution estimate x and a vector with the history of
residual norms. (This function is for demo only, not practical use.)
"""
function gmres(A,b,m)
    n = length(b)
    Q = zeros(n,m+1)
    Q[:,1] = b/norm(b)
    H = zeros(m+1,m)

    # Initial solution is zero.
    x = 0
    residual = [norm(b);zeros(m)]
    
    for j = 1:m
      # Next step of Arnoldi iteration.
      v = A*Q[:,j]
      for i = 1:j
          H[i,j] = dot(Q[:,i],v)
          v -= H[i,j]*Q[:,i]
      end
      H[j+1,j] = norm(v)
      Q[:,j+1] = v/H[j+1,j]

      # Solve the minimum residual problem.
      r = [norm(b); zeros(j)]
      z = H[1:j+1,1:j] \ r
      x = Q[:,1:j]*z
      residual[j+1] = norm( A*x - b )
    end

    return x,residual
end
```
````

## Convergence and restarting

Thanks to {prf:ref}`theorem-krylovmult`, minimization of $\|\mathbf{b}-\mathbf{A}\mathbf{x}\|$ over $\mathcal{K}_{m+1}$ includes minimization over $\mathcal{K}_m$. Hence the norm of the residual $\mathbf{r}_m = \mathbf{b} - \mathbf{A}\mathbf{x}_m$ (being the minimized quantity) cannot increase as the iteration unfolds.

```{index} convergence rate!linear
```

Unfortunately, making other conclusive statements about the convergence of GMRES is neither easy nor simple. [Our earlier demo](demos/gmres-intro.ipynb) shows the cleanest behavior: essentially linear convergence down to the range of machine epsilon. But it is possible for the convergence to go through phases of sublinear and superlinear convergence as well. There is a strong dependence on the spectrum of the matrix, a fact we state with more precision and detail in the next section.

```{index} GMRES!restarting
```

One of the practical challenges in GMRES is that as the dimension of the Krylov subspace grows, the number of new entries to be found in $\mathbf{H}_m$, and the total number of columns in $\mathbf{Q}$, also grow. Thus both the work and the storage requirements are quadratic in $m$, which can become intolerable in some applications. For this reason, GMRES is often used with {term}`restarting`.

Suppose $\hat{\mathbf{x}}$ is an approximate solution of $\mathbf{A}\mathbf{x}=\mathbf{b}$. Then if we set $\mathbf{x}=\mathbf{u}+\hat{\mathbf{x}}$, we have $\mathbf{A}(\mathbf{u}+\hat{\mathbf{x}}) = \mathbf{b}$, or $\mathbf{A}\mathbf{u} = \mathbf{b} - \mathbf{A}\hat{\mathbf{x}}$. The conclusion is that if we get an approximate solution and compute its residual $\mathbf{r}=\mathbf{b} - \mathbf{A}\hat{\mathbf{x}}$, then we need only to solve $\mathbf{A}\mathbf{u} = \mathbf{r}$ in order to get a correction to $\hat{\mathbf{x}}$.[^relativerestart]

[^relativerestart]: The new problem needs to be solved for accuracy relative to $\|\mathbf{b}\|$, *not* relative to $\|\mathbf{r}\|$.

Restarting guarantees a fixed upper bound on the per-iteration cost of GMRES. However, this bound comes at a price. Even though restarting preserves progress made in previous iterations, the Krylov space information is discarded and the residual minimization process starts again over low-dimensional choices. That can significantly retard or even stagnate the convergence. 

::::{prf:example} Julia demo
:class: demo
:label: demos-gmres-restart
{doc}`demos/gmres-restart`
::::

If the restart takes place before GMRES has entered a rapidly converging phase, the restarted GMRES can converge a great deal more slowly. However, the later iterations for the restarted versions should be much faster than those of pure GMRES, making an apples-to-apples comparison difficult.

<!-- There are other ways to avoid the growth in computational effort as the GMRES/Arnoldi iteration proceeds. Three of the more popular variations are abbreviated CGS, BiCGSTAB, and QMR, and these are also implemented in MATLAB. We do not describe them in this book. -->

<!-- 

\begin{exercises}
	\input{krylov/exercises/GMRES}
\end{exercises} -->
