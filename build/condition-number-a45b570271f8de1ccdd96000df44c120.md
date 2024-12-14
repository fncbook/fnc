(section-linsys-condition-number)=
# Conditioning of linear systems

```{index} condition number; of linear system
```

We are ready to consider the conditioning of solving the square linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$. In this problem, the data are $\mathbf{A}$ and $\mathbf{b}$, and the solution is $\mathbf{x}$. Both data and result are multidimensional, so we will use norms to measure their magnitudes.

The motivation for the definition of relative condition number in Chapter 1 was to quantify the response of the result to perturbations of the data. For simplicity, we start by allowing perturbations to $\mathbf{b}$ only while $\mathbf{A}$ remains fixed. 

Let $\mathbf{A}\mathbf{x}=\mathbf{b}$ be perturbed to

```{math}
  \mathbf{A}(\mathbf{x}+\mathbf{h}) = \mathbf{b}+\mathbf{d}.
```

The condition number should be the relative change in the solution divided by relative change in the data,

```{math}
  \frac{\quad\dfrac{\| \mathbf{h} \| }{\| \mathbf{x} \| }\quad}{\dfrac{\| \mathbf{d} \| }{\| \mathbf{b} \|}} 
  = \frac{\| \mathbf{h} \|\;\| \mathbf{b} \| }{\| \mathbf{d} \|\; \| \mathbf{x} \| }.
```

We can bound $\| \mathbf{h} \|$ in terms of $\| \mathbf{d} \|$:

```{math}
\begin{split}
  \mathbf{A}\mathbf{x} +  \mathbf{A} \mathbf{h} &= \mathbf{b} + \mathbf{d}, \\
  \mathbf{A} \mathbf{h} &= \mathbf{d},\\
  \mathbf{h} &= \mathbf{A}^{-1} \mathbf{d},\\
  \| \mathbf{h} \| &\le \| \mathbf{A}^{-1}\| \,\| \mathbf{d} \|,
\end{split}
```

where we have applied $\mathbf{A}\mathbf{x}=\mathbf{b}$ and {eq}`normineq1`. 
Since also $\mathbf{b}=\mathbf{A}\mathbf{x}$ implies $\| \mathbf{b} \|\le
\| \mathbf{A} \|\, \| \mathbf{x} \|$, we derive

```{math}
   \frac{\| \mathbf{h} \|\; \| \mathbf{b} \|}{\| \mathbf{d} \|\; \| \mathbf{x} \|} 
   \le \frac{\bigl(\| \mathbf{A}^{-1} \|\, \| \mathbf{d} \|\bigr)
    \bigl(\| \mathbf{A} \|\,\| \mathbf{x} \|\bigr)}{\| \mathbf{d} \|\,\| \mathbf{x} \|} 
    = \| \mathbf{A}^{-1}\| \, \| \mathbf{A} \|.
```

```{index} condition number; of a matrix
```

It is possible to show that this bound is tight, in the sense that the inequalities are in fact equalities for some choices of $\mathbf{b}$ and $\mathbf{d}$. This result motivates a new definition.

::::{prf:definition} Matrix condition number
The **matrix condition number** of an invertible square matrix $\mathbf{A}$ is 

```{math}
:label: mxcond
\kappa(\mathbf{A}) = \| \mathbf{A}^{-1}\| \, \| \mathbf{A} \|.
```

This value depends on the choice of norm; a subscript on $\kappa$ such as $1$, $2$, or $\infty$ is used if clarification is needed. If $\mathbf{A}$ is singular, we define $\kappa(\mathbf{A}) = \infty$.
::::

## Main result

The matrix condition number {eq}`mxcond` is equal to the condition number of solving a linear system of equations. Although we derived this fact above only for perturbations of $\mathbf{b}$, a similar statement holds when $\mathbf{A}$ is perturbed.

Using a traditional $\Delta$ notation for the perturbation in a quantity, we can write the following.

````{prf:theorem} Conditioning of linear systems
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

A condition number of 1 is the best we can hope for—in that case, the relative perturbation of the solution has the same size as that of the data.  A condition number of size $10^t$ indicates that in floating-point arithmetic, roughly $t$ digits are lost (i.e., become incorrect) in computing the solution $\mathbf{x}$. And if $\kappa(\mathbf{A}) > \epsilon_\text{mach}^{-1}$, then for computational purposes the matrix is effectively singular. 

(demo-condition-bound)=
```{prf:example}
```





```{index} ! Julia; cond
```
::::{grid} 1 1 2 2

:::{grid-item}
:columns: 7


Julia has a function `cond` to compute matrix condition numbers. By default, the 2-norm is used. As an example, the family of *Hilbert matrices* is famously badly conditioned. Here is the $6\times 6$  case.


:::
:::{card}
:columns: 5


Type `\kappa` followed by <kbd>Tab</kbd> to get the Greek letter $\kappa$.

:::
::::

```{code-cell}
A = [ 1/(i+j) for i in 1:6, j in 1:6 ]
κ = cond(A)
```

Because $\kappa\approx 10^8$, it's possible to lose nearly 8 digits of accuracy in the process of passing from $\mathbf{A}$ and $\mathbf{b}$ to $\mathbf{x}$. That fact is independent of the algorithm; it's inevitable once the data are expressed in finite precision. 

Let's engineer a linear system problem to observe the effect of a perturbation. We will make sure we know the exact answer.

```{code-cell}
x = 1:6
b = A*x
```

Now we perturb the system matrix and vector randomly by $10^{-10}$ in norm.

```{code-cell} 
ΔA = randn(size(A));  ΔA = 1e-10*(ΔA/opnorm(ΔA));
Δb = randn(size(b));  Δb = 1e-10*normalize(Δb);
```

We solve the perturbed problem using pivoted LU and see how the solution was changed.

```{code-cell}
new_x = ((A+ΔA) \ (b+Δb))
Δx = new_x - x
```

Here is the relative error in the solution.

```{code-cell}
@show relative_error = norm(Δx) / norm(x);
```

And here are upper bounds predicted using the condition number of the original matrix.

```{code-cell}
println("Upper bound due to b: $(κ*norm(Δb)/norm(b))")
println("Upper bound due to A: $(κ*opnorm(ΔA)/opnorm(A))")
```

Even if we didn't make any manual perturbations to the data, machine roundoff does so at the relative level of $\macheps$.

```{code-cell}
Δx = A\b - x
@show relative_error = norm(Δx) / norm(x);
@show rounding_bound = κ*eps();
```

Larger Hilbert matrices are even more poorly conditioned:

```{code-cell}
A = [ 1/(i+j) for i=1:14, j=1:14 ];
κ = cond(A)
```

Note that $\kappa$ exceeds $1/\macheps$. In principle we therefore may end up with an answer that has relative error greater than 100%.

```{code-cell}
rounding_bound = κ*eps()
```

Let's put that prediction to the test.

```{code-cell}
x = 1:14
b = A*x  
Δx = A\b - x
@show relative_error = norm(Δx) / norm(x);
```

As anticipated, the solution has zero accurate digits in the 2-norm.




## Residual and backward error

Suppose that $\mathbf{A}\mathbf{x}=\mathbf{b}$ and $\tilde{\mathbf{x}}$ is a computed estimate of the solution $\mathbf{x}$. The most natural quantity to study is the error, $\mathbf{x}-\tilde{\mathbf{x}}$. Normally we can't compute it because we don't know the exact solution. However, we can compute something related.

```{index} ! residual; of a linear system
```

::::{prf:definition} Residual of a linear system
For the problem $\mathbf{A}\mathbf{x}=\mathbf{b}$, the **residual** at a solution estimate $\tilde{\mathbf{x}}$ is 

```{math}
:label: residual
  \mathbf{r} = \mathbf{b} - \mathbf{A}\tilde{\mathbf{x}}.
```
::::

```{index} backward error; in a linear system
```

Obviously, a zero residual means that $\tilde{\mathbf{x}}=\mathbf{x}$, and we have the exact solution. What happens more generally? Note that $\mathbf{A}\tilde{\mathbf{x}}=\mathbf{b}-\mathbf{r}$. That is, $\tilde{\mathbf{x}}$ solves the linear system problem for a right-hand side that is changed by $-\mathbf{r}$. This is precisely what is meant by backward error.

Hence residual and backward error are the same thing for a linear system. What is the connection to the (forward) error? We can reconnect with {eq}`linsyscondb` by the definition $\mathbf{h} = \tilde{\mathbf{x}}-\mathbf{x}$, in which case 

$$\mathbf{d} = \mathbf{A}(\mathbf{x}+\mathbf{h})-\mathbf{b}=\mathbf{A}\mathbf{h} = -\mathbf{r}.$$ 

Thus {eq}`linsyscondb` is equivalent to

```{math}
:label: residualcond
  \frac{\| \mathbf{x}-\tilde{\mathbf{x}} \|}{\| \mathbf{x} \|} \le
  \kappa(\mathbf{A}) \frac{\| \mathbf{r} \|}{\| \mathbf{b} \|}.
```

Equation {eq}`residualcond` says that the gap between relative error and the relative residual is a multiplication by the matrix condition number. 

```{prf:observation}
When solving a linear system, all that can be expected is that the backward error, not the error, is small.
```

## Exercises

1. ⌨ Refer to {numref}`Demo {number} <demo-condition-bound>` for the definition of a Hilbert matrix. Make a table of the values of $\kappa(\mathbf{H}_n)$ in the 2-norm for $n=2,3,\ldots,16$. Speculate as to why the growth of $\kappa$ appears to slow down at $n=13$.

2. ⌨ The purpose of this problem is to verify, like in {numref}`Demo %s <demo-condition-bound>`, the error bound

    ```{math}
    \frac{\| \mathbf{x}-\tilde{\mathbf{x} \|}}{\| \mathbf{x} \|} \le \kappa(\mathbf{A})
    \frac{\| \mathbf{h} \|}{\| \mathbf{b} \|}.
    ```

    Here $\tilde{\mathbf{x}}$ is a numerical approximation to the exact solution $\mathbf{x}$, and $\mathbf{h}$ is an unknown perturbation caused by machine roundoff. We will assume that $\| \mathbf{d} \|/\| \mathbf{b} \|$ is roughly `eps()`.

    For each $n=10,20,\ldots,70$ let `A = matrixdepot("prolate",n,0.4)` and let $\mathbf{x}$ have components $x_k=k/n$ for $k=1,\ldots,n$. Define `b=A*x` and let $\tilde{\mathbf{x}}$ be the solution produced numerically by backslash.

    Make a table including columns for $n$, the condition number of $\mathbf{A}$, the observed relative error in $\tilde{\mathbf{x}}$, and the right-hand side of the inequality above. You should find that the inequality holds in every case.

3. ⌨ [Exercise 2.3.7](problem-ls-triangillcond) suggests that the solutions of linear systems

    ```{math}
    \mathbf{A} = \begin{bmatrix} 1 & -1 & 0 & \alpha-\beta & \beta \\ 0 & 1 & -1 &
      0 & 0 \\ 0 & 0 & 1 & -1 & 0 \\ 0 & 0 & 0 & 1 & -1  \\ 0 & 0 & 0 & 0 & 1
    \end{bmatrix}, \quad
    \mathbf{b} = \begin{bmatrix} \alpha \\ 0 \\ 0 \\ 0 \\ 1 \end{bmatrix}
    ```

    become less accurate as $\beta$ increases. Using $\alpha=0.1$ and $\beta=10,100,\ldots,10^{12}$, make a table with columns for $\beta$, $|x_1-1|$, and the condition number of the matrix.

4. ⌨ Let $\mathbf{A}_n$ denote the $(n+1)\times(n+1)$ version of the Vandermonde matrix in Equation {eq}`vandersystem` based on the equally spaced interpolation nodes $t_i=i/n$ for $i=0,\ldots,n$. Using the 1-norm, graph $\kappa(\mathbf{A}_n)$ as a function of $n$ for $n=4,5,6,\ldots,20$, using a log scale on the $y$-axis. (The graph is nearly a straight line.)
 

5. ⌨ The matrix $\mathbf{A}$ in {eq}`plu-stab-A` has unpivoted LU factors given in {eq}`plu-stab-LU` as a function of parameter $\epsilon$. For $\epsilon = 10^{-2},10^{-4},\ldots,10^{-10}$, make a table with columns for $\epsilon$, $\kappa(\mathbf{A})$, $\kappa(\mathbf{L})$, and $\kappa(\mathbf{U})$. (This shows that solution via unpivoted LU factorization is arbitrarily unstable.)
   
6. ✍  Define $\mathbf{A}_n$ as the $n\times n$ matrix $\displaystyle\begin{bmatrix}
      1 & -2 & & &\\
      & 1 & -2 & & \\
      & & \ddots & \ddots & \\
      & & & 1 & -2 \\
      & & & & 1
    \end{bmatrix}.$

    **(a)** Write out $\mathbf{A}_2^{-1}$ and $\mathbf{A}_3^{-1}$.
  
    **(b)** Write out $\mathbf{A}_n^{-1}$ in the general case $n>1$. (If necessary, look at a few more cases in Julia until you are certain of the pattern.) Make a clear argument why it is correct.

    **(c)** Using the $\infty$-norm, find $\kappa(\mathbf{A}_n)$.

7. ✍ **(a)** Prove that for $n\times n$ nonsingular matrices $\mathbf{A}$ and $\mathbf{B}$, $\kappa(\mathbf{A}\mathbf{B})\le \kappa(\mathbf{A})\kappa(\mathbf{B})$.

    **(b)** Show by means of an example that the result of part (a) cannot be an equality in general.

8. ✍  Let $\mathbf{D}$ be a diagonal $n\times n$ matrix, not necessarily invertible. Prove that in the 1-norm,

    ```{math}
    \kappa(\mathbf{D}) = \frac{\max_i |D_{ii}|}{\min_i |D_{ii}|}.
    ```

    (Hint: See [Exercise 2.7.10](problem-norms-diagnorm).)

