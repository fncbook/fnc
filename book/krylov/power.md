# Power iteration

```{index} sparse matrix
```

Given that matrix-vector multiplication is fast for sparse matrices, let's see what we might accomplish with only that at our disposal.

::::{sidebar} Demo
:class: demo
{doc}`demos/power-one`
::::

The [results of the demo](demos/power-one) lead to some speculation. If it holds exactly, the equation $\mathbf{A}\mathbf{x} = \mathbf{x}$ is a special case of the eigenvalue condition $\mathbf{A}\mathbf{x} = \lambda \mathbf{x}$. It appears, then, that as $k\to\infty$, $\mathbf{A}^k \mathbf{x}$ converges to an eigenvector of $\mathbf{A}$ for the eigenvalue $\lambda=1$.

## Dominant eigenvalue

```{index} eigenvalue; dominant
```
Analysis of matrix powers is fairly straightforward, if the matrix is diagonalizable. Let $\mathbf{A}$ be any diagonalizable $n\times n$ matrix having eigenvalues $\lambda_1,\ldots,\lambda_n$ and corresponding linearly independent eigenvectors $\mathbf{v}_1,\ldots,\mathbf{v}_n$. Furthermore, suppose the eigenvalues (possibly after renumbering) are such that

:::{math}
:label: evorder
|\lambda_1| > |\lambda_2| \ge |\lambda_3| \ge \cdots \ge |\lambda_n|.
:::

Given {eq}`evorder` we say that $\lambda_1$ is the {term}`dominant eigenvalue`. This was the case with $\lambda_1=1$ for $\mathbf{A}$ in [the demo above](demos/power-one).

Now let the initially chosen $\mathbf{x}$ be expressed as a linear combination in the eigenvector basis:

:::{math}
\mathbf{x} = c_1 \mathbf{v}_{1} + c_2 \mathbf{v}_{2} + \cdots + c_n\mathbf{v}_{n}.
:::

Observe that

\begin{align*}
\mathbf{A}\mathbf{x} &= c_1 \mathbf{A}\mathbf{v}_{1} + c_2 \mathbf{A}\mathbf{v}_{2} + \cdots + c_n \mathbf{A}\mathbf{v}_{n}\\
&= c_1 \lambda_1 \mathbf{v}_{1} + c_2 \lambda_2  \mathbf{v}_{2} +
\cdots + c_n \lambda_n \mathbf{v}_{n}.
\end{align*}

If we apply $\mathbf{A}$ repeatedly, then

::::{math}
:label: powerAkx0
\begin{split}
  \mathbf{A}^k\mathbf{x} &= \lambda_1^k c_1 \mathbf{v}_{1} + \lambda_2^k c_2
			\mathbf{v}_{2} + \cdots + \lambda_n^k c_n \mathbf{v}_{n}  \\
	&= \lambda_1^k \left[ c_1 \mathbf{v}_{1} +
		\left(\frac{\lambda_2}{\lambda_1}\right) ^k c_2
		\mathbf{v}_{2} + \cdots + \left(\frac{\lambda_n}{\lambda_1}\right)^k
		c_n \mathbf{v}_{n} \right].
\end{split}
::::

Since $\lambda_1$ is dominant, we conclude that in any norm,

:::{math}
:label: poweriterconverge
\left\| \frac{ \mathbf{A}^k\mathbf{x}}{\lambda_1^k}
- c_1\mathbf{v}_1\right\| \le |c_2|\cdot\left|\frac{\lambda_2}{\lambda_1}\right| ^k
\| \mathbf{v}_{2} \| + \cdots +  |c_n|\cdot\left|\frac{\lambda_n}{\lambda_1}\right|^k
\| \mathbf{v}_{n} \| \rightarrow 0 \text{ as $k\rightarrow \infty$}.
:::

As long as $c_1\neq 0$, $\mathbf{A}^k\mathbf{x}$ eventually is almost parallel to the dominant eigenvector.[^zeromeasure]

[^zeromeasure]: If $\mathbf{x}$ is chosen randomly, the odds that $c_1=0$ are mathematically zero.

For algorithmic purposes, it is important to interpret $\mathbf{A}^k\mathbf{x}$ as $\mathbf{A}\bigl( \cdots\bigl( \mathbf{A} (\mathbf{A}\mathbf{x})\bigl) \cdots\bigl)$, i.e. as repeated applications of $\mathbf{A}$ to a vector. This interpretation allows us to fully exploit any sparsity of $\mathbf{A}$, something which is [not preserved](demos/structure-fill.ipynb) by taking a matrix power $\mathbf{A}^k$ explicitly.

## Power iteration

```{index} power iteration
```
An important technicality separates us from an algorithm: unless $|\lambda_1|=1$, the factor $\lambda_1^k$ tends to make $\|\mathbf{A}^k\mathbf{x}\|$ either very large or very small. To make a practical algorithm, we alternate matrix-vector multiplication with a renormalization of the vector. This algorithm is known as the {term}`power iteration`.

1. Choose $\mathbf{x}_1$.
2. For $k=1,2,\ldots$
	::::{math}
    :label: poweriter
    \begin{split}
    \mathbf{y}_k &= \mathbf{A} \mathbf{x}_k, \\
    \alpha_k &= \frac{1}{y_{k,m}}, \text{ where } |y_{k,m}|=\|{\mathbf{y}_k} \|_\infty,\\
    \mathbf{x}_{k+1} &= \alpha_k \mathbf{y}_k.   
    \end{split}
	::::

Our notation here uses $y_{k,m}$ to mean the $m$th component of $\mathbf{y}_k$, or $\mathbf{e}_m^T\mathbf{y}_k$. Note that now $\| \mathbf{x}_{k+1}\|_\infty=1$. We can write

:::{math}
:label: powernorm
\mathbf{x}_{k} = (\alpha_1 \alpha_2 \cdots \alpha_k ) \mathbf{A}^k \mathbf{x}_{1}.
:::

Thus the renormalization step modifies {eq}`powerAkx0` and {eq}`poweriterconverge` only slightly.

So far we have discussed only eigenvector estimation. However, if $\mathbf{x}_k$ is nearly a dominant eigenvector of $\mathbf{A}$, then $\mathbf{A}\mathbf{x}_k$ is nearly $\lambda_1\mathbf{x}_k$, and we can take the ratio $\gamma_k=y_{k,m}/x_{k,m}$ as an eigenvalue estimate. In fact, revisiting {eq}`powerAkx0`, the extra $\alpha_j$ normalization factors cancel in the ratio, and, after some simplification, we get

:::{math}
:label: poweriterratio
\gamma_k = \frac{y_{k,m}}{x_{k,m}} = \lambda_1
\frac{1+r_2^{k+1} b_2 + \cdots +  r_n^{k+1} b_n}{1+r_2^{k} b_2 +  \cdots +  r_n^{k} b_n},
:::

where $r_j=\lambda_j/\lambda_1$ and the $b_j$ are constants. By assumption {eq}`evorder`, each $r_j$ satisfies $|r_j|<1$, so we see that $\gamma_k\rightarrow \lambda_1$ as $k\rightarrow\infty$.

Here is our implementation of power iteration.

(function-poweriter)=
````{proof:function} poweriter
**Power iteration to find a dominant eigenvalue.**

```{code-block} julia
:lineno-start: 1
"""
poweriter(A,numiter)

Perform `numiter` power iterations with the matrix `A`, starting
from a random vector, and return a vector of eigenvalue estimates
and the final eigenvector approximation.
"""
function poweriter(A,numiter)
    n = size(A,1)
    x = normalize(randn(n),Inf)
    gamma = zeros(numiter)
    for k = 1:numiter
      y = A*x
      normy,m = findmax(abs.(y))
      gamma[k] = y[m]/x[m]
      x = y/y[m]
    end

    return gamma,x
end
```
````

The {term}`findmax` function used above returns both the largest element of a vector and its location in the vector. Observe that the only use of $\mathbf{A}$ is to find the matrix-vector product $\mathbf{A}\mathbf{x}$, which makes exploitation of the sparsity of $\mathbf{A}$ automatic.

## Convergence

Let's examine the terms in the numerator and denominator of {eq}`poweriterratio` more carefully:

:::{math}
:label: powerleftover
\begin{split}
r_2^{k} b_2 +  \cdots +  r_n^{k} b_n &= r_2^k \left[ b_2 + \left( \frac{r_3}{r_2} \right)^kb_3 + \cdots + \left( \frac{r_n}{r_2} \right)^kb_n \right] \\
&= r_2^k \left[ b_2 + \left( \frac{\lambda_3}{\lambda_2} \right)^kb_3 + \cdots + \left( \frac{\lambda_n}{\lambda_2} \right)^kb_n \right].
\end{split}
:::

At this point we'll introduce an extra assumption,

:::{math}
:label: evorder2
|\lambda_2| > |\lambda_3| \ge \cdots \ge |\lambda_n|.
:::

This condition isn't strictly necessary, but it simplifies the following statements considerably, because now it's clear that the quantity in {eq}`powerleftover` approaches $b_2 r_2^k$ as $k\rightarrow \infty$.

Next we estimate {eq}`poweriterratio` for large $k$, using a geometric series expansion for the denominator to get

:::{math}
:label: powerest
\begin{split}
\gamma_k & \to \lambda_1 \left( 1+b_2 r_2^{k+1} \right) \left( 1 - b_2 r_2^{k} + O(r_2^{2k}) \right), \\
\gamma_k - \lambda_1 &\to \lambda_1 b_2 (  r_2 - 1 ) r_2^{k}.
\end{split}
:::

```{index} convergence rate!linear
```
This is linear convergence with factor $r_2$. That is,

:::{math}
:label: poweriterconv
\frac{\gamma_{k+1} - \lambda_1}{\gamma_{k}-\lambda_1} \rightarrow r_2 = \frac{\lambda_2}{\lambda_1} \quad \text{as } k\rightarrow \infty.
:::

The error in the eigenvalue estimates $\gamma_k$ of power iteration is reduced asymptotically by a constant factor $\lambda_2/\lambda_1$ on each iteration.

::::{sidebar} Demo
:class: demo
{doc}`demos/power-iter`
::::

The practical utility of {eq}`poweriterconv` is limited, because if we knew $\lambda_1$ and $\lambda_2$, we wouldn't be running the power iteration in the first place! Sometimes it's possible to find estimates or bounds of the ratio. But for the most part we just find it a valuable theoretical statement of how power iteration should converge.

<!-- \begin{exercises}
    \input{krylov/exercises/PowerIteration}
\end{exercises} -->

