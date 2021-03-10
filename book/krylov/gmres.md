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

Compare the graph in {prf:ref}`demos-gmres-intro`  to the one in {prf:ref}`demos-subspace-unstable`. Both start with the same linear convergence, but only the version using Arnoldi avoids the instability created by the poor Krylov basis.

A basic implementation of GMRES is given below.

(function-gmres)=
````{proof:function} gmres
**GMRES for a linear system**

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

```{index} convergence rate; linear
```

Unfortunately, making other conclusive statements about the convergence of GMRES is neither easy nor simple. {prf:ref}`demos-gmres-intro` shows the cleanest behavior: essentially linear convergence down to the range of machine epsilon. But it is possible for the convergence to go through phases of sublinear and superlinear convergence as well. There is a strong dependence on the spectrum of the matrix, a fact we state with more precision and detail in the next section.

```{index} GMRES; restarting
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

There are other ways to avoid the growth in computational effort as the GMRES/Arnoldi iteration proceeds. Three of the more popular variations are abbreviated CGS, BiCGSTAB, and QMR. We do not describe them in this book.

## Exercises

1. ✍ (See also [an earlier exercise](problem-krylovpermute).) Consider the linear system with
  
    $$
    \mathbf{A}=\displaystyle 
    \begin{bmatrix}
      0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \\ 1 & 0 & 0 & 0
    \end{bmatrix}, \qquad \mathbf{b}=\mathbf{e}_1.
    $$

    **(a)** Find the exact solution by inspection.

    **(b)** Find the GMRES approximate solutions $\mathbf{x}_m$ for $m=1,2,3,4$. 


    ::::{only} solutions
    The Krylov vectors are $b=e_1$, $Ab=e_4$, $A^2b=e_3$,
      $A^3b=e_2$. 
    ::::

2. ✍ (Continuation of [an earlier exercise](problem-matrixpolykrylov).) Show that if $\mathbf{x}_m\in\mathcal{K}_m$, then the residual $\mathbf{b}-\mathbf{A}\mathbf{x}_m$ is equal to $q(\mathbf{A})\mathbf{b}$, where $q$ is a polynomial of degree at most $m$ and $q(0)=1$. (This fact is a key one for many convergence results.) 

3. ✍ Explain why GMRES, in exact arithmetic, converges to the true solution in $n$ iterations for an $n\times n$ matrix if $\operatorname{rank}(\mathbf{K}_n)=n$. (Hint: Consider how the algorithm is defined from first principles.) 

4. ⌨ Let $\mathbf{A}$ be the $n\times n$ tridiagonal matrix
  
    \begin{bmatrix}
      -4 & 1      &        &        &   \\
      1  & -4     & 1      &        &   \\
         & \ddots & \ddots & \ddots &   \\
         &        & 1      & -4     & 1 \\
         &        &        & 1      & -4 
    \end{bmatrix}
  
    and let the $n$-vector $\mathbf{b}$ have entries $b_i=i/n$. For $n=8,16,32,64$, run {numref}`Function {number}<function-gmres>` for $m=n/2$ iterations. On one graph plot $\|\mathbf{r}_k\|/\|\mathbf{b}\|$ for all the cases. How does the convergence rate of GMRES seem to depend on $n$?  

    (problem-gmres-surround)=
5. ⌨  In this problem you will see the strong effect the eigenvalues of the matrix may have on GMRES convergence. Let 
   
    $$
    \mathbf{B}=
    \begin{bmatrix}
      1 & & & \\
      & 2 & & \\
      & & \ddots & \\
      & & & 100
    \end{bmatrix},
    $$ 
    
    let $\mathbf{I}$ be a $100\times 100$ identity, and let $\mathbf{Z}$ be a $100\times 100$ matrix of zeros. Also let $\mathbf{b}$ be a $200\times 1$ vector of ones. 

    **(a)** Let $\mathbf{A} = \begin{bmatrix}
        \mathbf{B} & \mathbf{I} \\ \mathbf{Z} & \mathbf{B}
      \end{bmatrix}.$ What are its eigenvalues (no computer required here)? Apply `gmres` with tolerance $10^{-10}$ for 100 iterations without restarts, and plot the residual convergence. 
    
    **(b)** Repeat part (a) with restarts every 20 iterations. 
    
    **(c)** Now let $\mathbf{A} = \begin{bmatrix}
        \mathbf{B} & \mathbf{I} \\ \mathbf{Z} & -\mathbf{B}
      \end{bmatrix}.$ What are its eigenvalues? Repeat part~(a). Which matrix is more difficult for GMRES? 

      ::::{only} solutions
      %% (a)
      B = diag(1:100);
      A = [B eye(100); zeros(100) B];
      b = ones(200,1);
      [~,~,~,~,rv] = gmres(A,b,100,1e-10,1);
      semilogy(rv)
      %% (b) 
      [~,~,~,~,rv] = gmres(A,b,20,1e-10,5);
      semilogy(rv)
      %% (c)
      A = [B eye(100); 0*I -B];
      [~,~,~,~,rv] = gmres(A,b,100,1e-10,1);
      semilogy(rv)    % much less convergence  
      ::::

6. ⌨ (Continuation of [an earlier exercise](problem-power-lumpmembraneinveig).) We again consider the $n^2\times n^2$ sparse matrix defined by `FNC.poisson(n)`.
The solution of $\mathbf{A}\mathbf{x}=\mathbf{b}$ may be interpreted as the deflection of a lumped membrane in response to a constant load represented by $\mathbf{b}$.
    
    **(a)** For  $n=10,15,20,25$, let $\mathbf{b}$ be the vector of $n^2$ ones and apply {numref}`Function {number} <function-gmres>` for 50 iterations. On one semi-log graph, plot the four convergence curves $\|\mathbf{r}_m\|/\|\mathbf{b}\|$.

    **(b)** For the case $n=25$ use "plot(1:n,1:n,reshape(x,25,25))" to plot the solution, which should look physically plausible. 

    ::::{only} solutions
    A = @(n) n^2*gallery('poisson',n);
    clf
    nValues = 10:5:25;
    for nI = 1:length(nValues)
        n = nValues(nI)
        b = ones(n^2,1);
        [x,r] = arngmres(A(n),b,50);
        semilogy(0:50,r/norm(b),'-o'), hold on
    end
    %%
    clf
    mesh(reshape(x,25,25))
    ::::

    (problem-helmhotzmatrix)=
7. ⌨  Consider a matrix that arises from the Helmholtz equation for wave propagation. The matrix can be specified using 

    ```julia
    A = FNC.poisson(n) - k^2*I;
    ```
  
    Let $n=20$ and use the two parameter values $k=1.3\pi$ and $k=1.5\pi$. Use `eigs` to find the five largest and five smallest eigenvalues in each case. (See {prf:ref}`demos-structure-linalg` for examples of its usage.) What happens to the spectrum of the matrix?  Specifically, does it remain definite?

8. ✍ Recall {eq}`gmresproblem`, defining the stable form of the GMRES minimization problem.

    **(a)** Show that $\mathbf{b}=\mathbf{Q}_{m+1}\tilde{\mathbf{b}}$ for an $(m+1)\times 1$ vector $\tilde{\mathbf{b}}$ that is nonzero only in its first component.

    **(b)** Show that the quantity $\| \mathbf{A}\mathbf{Q}_m\mathbf{z}-\mathbf{b} \|$ is equal to $\| \mathbf{H}_m\mathbf{z} - \tilde{\mathbf{b} \|}$.
      
    **(c)** Modify {numref}`Function {number}<function-gmres>`, replacing line 33 with the solution of a linear least-squares problem of size $(m+1)\times m$.


