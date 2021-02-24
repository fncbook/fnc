# Conditioning of linear systems

```{index} condition number; of linear system
```

We are ready to consider the conditioning of solving the square linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$. Recall that the condition number is the relative change in the solution divided by a relative change in the data. In this case the data are $\mathbf{A}$ and $\mathbf{b}$, and the solution is $\mathbf{x}$.

For simplicity we start by allowing perturbations to $\mathbf{b}$ only while $\mathbf{A}$ remains unchanged. Let $\mathbf{A}\mathbf{x}=\mathbf{b}$ be perturbed to

```{math}
  \mathbf{A}(\mathbf{x}+\mathbf{h}) = \mathbf{b}+\mathbf{d}.
```

We seek to bound $\| \mathbf{h} \|$ in terms of $\| \mathbf{d} \|$:

```{math}
\begin{split}
  \mathbf{A}\mathbf{x} +  \mathbf{A} \mathbf{h} &= \mathbf{b} + \mathbf{d} \\
  \mathbf{A} \mathbf{h} &= \mathbf{d}\\
  \mathbf{h} &= \mathbf{A}^{-1} \mathbf{d}\\
  \| \mathbf{h} \| &\le \| \mathbf{A}^{-1}\| \,\| \mathbf{d} \|,
\end{split}
```

where we have used $\mathbf{A}\mathbf{x}=\mathbf{b}$ and {eq}`normineq1`.
Since furthermore $\mathbf{b}=\mathbf{A}\mathbf{x}$ and therefore $\| \mathbf{b} \|\le
\| \mathbf{A} \|\, \| \mathbf{x} \|$, we derive

```{math}
  \frac{\quad\dfrac{\| \mathbf{h} \|}{\| \mathbf{x} \|}\quad}{\dfrac{\| \mathbf{d} \|}{\| \mathbf{b} \|}} = \frac{\| \mathbf{h} \|\;
    \| \mathbf{b} \|}{\| \mathbf{d} \|\; \| \mathbf{x} \|} \le
  \frac{\bigl(\| \mathbf{A}^{-1} \|\, \| \mathbf{d} \|\bigr)
    \bigl(\| \mathbf{A} \|\,\| \mathbf{x} \|\bigr)}{\| \mathbf{d} \|\,\| \mathbf{x} \|} =
  \| \mathbf{A}^{-1}\| \, \| \mathbf{A} \|.
```

```{index} condition number; of a matrix
```

It is possible to show that this bound is tight, in the sense that the inequalities are in fact equalities for some choices of $\mathbf{b}$ and $\mathbf{d}$. Motivated by the definition of the condition number as the ratio of relative changes in solution and data, we define the {term}`matrix condition number`

```{math}
:label: mxcond
\kappa(\mathbf{A}) = \| \mathbf{A}^{-1}\| \, \| \mathbf{A} \|.
```

Note that $\kappa(\mathbf{A})$ depends on the choice of norm; a subscript on $\kappa$ such
as $1$, $2$, or $\infty$ is used if clarification is needed. 

## Main result

```{margin}
The matrix condition number is equal to the condition number of solving a linear system of equations.
```

The matrix condition number {eq}`mxcond` is equal to the condition number of solving a linear system of equations. Although we derived this fact only for perturbations of $\mathbf{b}$, a similar statement holds when $\mathbf{A}$ is perturbed.

Using a traditional $\Delta$ notation for the perturbation in a quantity, we can write the following.

````{proof:observation}

If $\mathbf{A}(\mathbf{x} + \Delta \mathbf{x}) = \mathbf{b} + \Delta \mathbf{b}$, then

```{math}
:label: linsyscondb
\frac{\| \Delta \mathbf{x} \|}{\| \mathbf{x} \|} \le \kappa(\mathbf{A}) \frac{\| \Delta \mathbf{b} \|}{\| \mathbf{b} \|}.
```

If $(\mathbf{A}+\Delta \mathbf{A}) (\mathbf{x} + \Delta \mathbf{x}) = \mathbf{b}$, then

```{math}
:label: linsyscondA
\frac{\| \Delta \mathbf{x} \|}{\| \mathbf{x} \|} \le \kappa(\mathbf{A}) \frac{\| \Delta \mathbf{A} \|}{\| \mathbf{A} \|},
```

in the limit $\| \Delta \mathbf{A} \| \to 0$.
````

Note that for any induced matrix norm,

```{math}
  1 = \| \mathbf{I} \| = \| \mathbf{A} \mathbf{A}^{-1} \| \le \| \mathbf{A} \|\, \| \mathbf{A}^{-1} \| = \kappa(\mathbf{A}).
```

A condition number of 1 is the best we can hope for—in that case, the relative perturbation of the solution has the same size as that of the data.  A condition number of size $10^t$ indicates that in floating point arithmetic, roughly $t$ digits are lost (i.e., become incorrect) in computing the solution $\mathbf{x}$.

```{margin}
If $\kappa(\mathbf{A}) > \epsilon_\text{mach}^{-1}$, then for computational purposes $\mathbf{A}$ is singular.
```

```{sidebar} Demo
:class: demo
{doc}`demos/condition-bound`
```

If $\kappa(\mathbf{A}) > \epsilon_\text{mach}^{-1}$, then for computational purposes the matrix is singular. If $\mathbf{A}$ is exactly singular, it is customary to say that $\kappa(\mathbf{A})=\infty$.

## Residual and backward error

```{index} residual
```

Suppose that $\mathbf{A}\mathbf{x}=\mathbf{b}$ and $\tilde{\mathbf{x}}$ is a computed estimate of the solution $\mathbf{x}$. The most natural quantity to study is the error, $\mathbf{x}-\tilde{\mathbf{x}}$. Normally we can't compute it because we don't know the exact solution. However, we can certainly compute the {term}`residual`, defined as

```{math}
  :label: residual
  \mathbf{r} = \mathbf{b} - \mathbf{A}\tilde{\mathbf{x}}.
```

```{index} backward error; in linear system
```

Obviously a zero residual means that $\tilde{\mathbf{x}}=\mathbf{x}$ and we get the exact solution. What does a "small" residual mean? Note that $\mathbf{A}\tilde{\mathbf{x}}=\mathbf{b}-\mathbf{r}$, so $\tilde{\mathbf{x}}$ solves the linear system problem for a right-hand side that is changed by $-\mathbf{r}$. This is precisely what is meant by backward error: the perturbation from the original problem to the one that is solved exactly.

But does a small residual mean that the error is also small? We can reconnect with {eq}`linsyscondb` by the definition $\mathbf{h} = \tilde{\mathbf{x}}-\mathbf{x}$, in which case $\mathbf{d} = \mathbf{A}(\mathbf{x}+\mathbf{h})-\mathbf{b}=\mathbf{A}\mathbf{h} = -\mathbf{r}$. Hence {eq}`linsyscondb` is equivalent to

```{math}
  :label: residualcond
  \frac{\| \mathbf{x}-\tilde{\mathbf{x} \|}}{\| \mathbf{x} \|} \le
  \kappa(\mathbf{A}) \frac{\| \mathbf{r} \|}{\| \mathbf{b} \|}.
```

```{margin}
When solving a linear system, all that can be expected is that the backward error, not the error, is small.
```

Equation {eq}`residualcond` says that the relative error can be much larger than the relative residual when the matrix condition number is large. To put it another way: When solving a linear system, all that can be expected is that the backward error, not the error, be small.

## Exercises

1. ⌨ A **Hilbert matrix** is a square matrix whose $(i,j)$ entry is $1/(i+j-1)$. The $n\times n$ version $\mathbf{H}_n$ can be generated in Julia using

    ``` julia
    H = [ 1/(i+j-1) for i in 1:n, j in 1:n ]
    ```

    The condition number of a Hilbert matrix grows very rapidly as a function of $n$, showing that even simple, small linear systems can be badly conditioned. 

    Make a table of the values of $\kappa(\mathbf{H}_n)$ in the 2-norm for $n=2,3,\ldots,16$. Why does the growth of $\kappa$ appear to slow down at $n=13$?

2. ⌨ The purpose of this problem is to verify, like in {doc}`demos/condition-bound`, the error bound

    ```{math}
    \frac{\| \mathbf{x}-\tilde{\mathbf{x} \|}}{\| \mathbf{x} \|} \le \kappa(\mathbf{A})
    \frac{\| \mathbf{h} \|}{\| \mathbf{b} \|}.
    ```

    Here $\tilde{\mathbf{x}}$ is a numerical approximation to the exact solution $\mathbf{x}$, and $\mathbf{h}$ is an unknown perturbation caused by
    machine roundoff. We will assume that $\| \mathbf{d} \|/\| \mathbf{b} \|$ is roughly of size `eps()`.

    For this problem you will need the `MatrixDepot` package, which can be installed and loaded via 

    ``` julia
    import Pkg; Pkg.add("MatrixDepot")
    using MatrixDepot
    ```

    For each $n=10,20,\ldots,70$ let `A = matrixdepot("prolate",n,0.4)` and let $\mathbf{x}$ have components $x_k=k/n$ for $k=1,\ldots,n$. Define `b=A*x` and let $\tilde{\mathbf{x}}$ be the solution produced numerically by backslash.

    Make a table including columns for $n$, the relative error in $\tilde{\mathbf{x}}$, the condition number of $\mathbf{A}$, and the right-hand side of the inequality above. You should find that the inequality holds in every case.

      %%
      % We assume that the relative change in $b$ is |eps|; that is,
      % $\|\Delta b\|/\|b\|=$|eps|.
      %%
      <!-- n = [10:10:70];
      relerr = []; condA = []
      for k = 1:length(n)
          A = gallery('prolate',n(k),0.4);
          x_exact = [1:n(k)]'/n(k);
          b = A*x_exact;
          x = A\b;
          relerr = [relerr norm(x-x_exact)/norm(x)];
          condA = [condA cond(A)];
      end -->

      %%
      % Output:
      <!-- disp(' ')    %  a blank line
      disp('  Problem 2.7.1                         ')
      disp(sprintf(' n     rel err        kappa(A)        rhs   '))
      disp('-----------------------------------------------')
      for k=1:length(n)
          disp(sprintf('%2.0f   %6.4e    %6.4e    %6.4e',...
          n(k),relerr(k),condA(k),condA(k)*eps))
      end -->

      %%
      % The bound works in every case as shown by the results below.

3. ⌨ An [earlier problem](problem-triangillcond) asked you to solve systems

    ```{math}
    \mathbf{A} = \begin{bmatrix} 1 & -1 & 0 & \alpha-\beta & \beta \\ 0 & 1 & -1 &
      0 & 0 \\ 0 & 0 & 1 & -1 & 0 \\ 0 & 0 & 0 & 1 & -1  \\ 0 & 0 & 0 & 0 & 1
    \end{bmatrix}, \quad
    \mathbf{b} = \begin{bmatrix} \alpha \\ 0 \\ 0 \\ 0 \\ 1 \end{bmatrix}
    ```

    with $\alpha=0.1$ and $\beta=10,100,\ldots,10^{12}$. Again make a table of $\beta$ and $|x_1-1|$, and add a column for the condition numbers of these matrices.

4. ⌨ The condition number of the formulation of polynomial interpolation suggested in {doc}`polyinterp` grows rapidly as the degree of the polynomial increases. Let $\mathbf{A}_n$ denote the $(n+1)\times(n+1)$ version of the Vandermonde matrix in equation {eq}`vandersystem` based on the equally spaced interpolation nodes $t_i=i/n$ for $i=0,\ldots,n$.
  
    **(a)** Using the 1-norm, graph $\kappa(\mathbf{A}_n)$ as a function of $n$ for $n=4,5,6,\ldots,20$, using a log scale on the $y$-axis. (The graph is nearly a straight line.)

    **(b)** Show that if $\log f(n)$ is a linear function of $n$ with positive slope, then $f(n) = C\alpha^n$ for constants $C>0$ and $\alpha>1$.
  
    <!-- n_ = 4:20;
    kappa_ = [];
    for n = n_
        t = linspace(0,1,n+1);
        A = vander(t);
        kappa_ = [kappa_;cond(A,1)];
    end
    semilogy(n_,kappa_,'o-') -->

5. ✍  Define $\mathbf{A}_n$ as the $n\times n$ matrix $\displaystyle
    \begin{bmatrix}
      1 & -2 & & &\\
      & 1 & -2 & & \\
      & & \ddots & \ddots & \\
      & & & 1 & -2 \\
      & & & & 1
    \end{bmatrix}.$

    **(a)** Write out $\mathbf{A}_2^{-1}$ and $\mathbf{A}_3^{-1}$.
  
    **(b)** Write out $\mathbf{A}_n^{-1}$ in the general case $n>1$. (If necessary, look at a few more cases in Julia until you are certain of the pattern.) Make a clear argument why it is correct.

    **(c)** Using the $\infty$-norm, find $\kappa(\mathbf{A}_n)$.

6. ✍ **(a)** Prove that for $n\times n$ nonsingular matrices $\mathbf{A}$ and $\mathbf{B}$, $\kappa(\mathbf{A}\mathbf{B})\le \kappa(\mathbf{A})\kappa(\mathbf{B})$.

    **(b)** Show by means of an example that the result of part (a) cannot be an equality in general.

    %% (a)
    %
    % $$\kappa(AB) = \|AB\|\cdot\|(AB)^{-1}\| \le \|A\|\cdot\|B\|\cdot\|B^{-1}A^{-1}\|\le \kappa(A)\kappa(B).$$
    %% (b)
    % Virtually all (nonsingular) matrix pairs will have a true inequality. For example,
    <!-- A = rand(5);  B = magic(5);
    cond(A*B) - cond(A)*cond(B) -->

7. ✍  Let $\mathbf{D}$ be a diagonal $n\times n$ matrix, not necessarily invertible. Prove that in the 2-norm,

    ```{math}
    \kappa(\mathbf{D}) = \frac{\max_i |D_{ii}|}{\min_i |D_{ii}|}.
    ```

    (Hint: See [this previous problem](problem-diagnorm).)

    % The determinant of $D$ is the product of the $d_{ii}$, which is zero if and only if one of the diagonal entries is zero. Say $D$ is singular. Then $min_i |d_{ii}|=0$, so $\kappa(D)=\infty$, as required.
    %%
    % If $D$ is invertible, then $D^{-1}={\rm diag}\,(d_{11}^{-1},\ldots,d_{nn}^{-1})$ and $\|D^{-1}\|_2=\max |d_{ii}^{-1}| = \left[ \min |d_{ii}| \right]^{-1}$, which leads directly to the result.
