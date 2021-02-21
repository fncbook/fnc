# Dimension reduction

```{index} singular value decomposition
```
```{index} dimension reduction
```

The SVD has another important property that proves very useful in a variety of applications. Let $\mathbf{A}$ be a real $m\times n$ matrix with SVD $\mathbf{A}=\mathbf{U}\mathbf{S}\mathbf{V}^T$ and (momentarily) $m\ge n$. Another way of writing the thin form of the SVD is

::::{math}
:label: svdsum
\begin{split}
  \mathbf{A} = \widehat{\mathbf{U}}\widehat{\mathbf{S}}\mathbf{V}^T &=
  \begin{bmatrix}
    \mathbf{u}_1 & \mathbf{u}_2 & \cdots & \mathbf{u}_n
  \end{bmatrix}
  \begin{bmatrix}
    \sigma_1 & & \\
    & \ddots & \\
    & & \sigma_n
  \end{bmatrix}
        \begin{bmatrix}
          \mathbf{v}_1^T \\ \vdots \\ \mathbf{v}_n^T
        \end{bmatrix}  \notag \\
  &=
  \begin{bmatrix}
    \sigma_1\mathbf{u}_1  & \cdots & \sigma_n\mathbf{u}_n
  \end{bmatrix}
                               \begin{bmatrix}
                                  \mathbf{v}_1^T \\ \vdots \\ \mathbf{v}_n^T
                               \end{bmatrix} \notag \\
  &= \sigma_1 \mathbf{u}_{1}\mathbf{v}_{1}^T + \cdots + \sigma_r \mathbf{u}_{r}\mathbf{v}_{r}^T = \sum_{i=1}^r \sigma_i \mathbf{u}_{i}\mathbf{v}_{i}^T,
\end{split}
::::

where $r$ is the rank of $\mathbf{A}$. The final formula also holds for the case $m<n$.

Each outer product $\mathbf{u}_{i}\mathbf{v}_{i}^T$ is a rank-1 matrix of unit norm. Thanks to the ordering of singular values, then, equation {eq}`svdsum` expresses $\mathbf{A}$ as a sum of decreasingly important contributions. This motivates the definition, for $1\le k \le r$,

```{math}
:label: svdlowrank
\mathbf{A}_k = \sum_{i=1}^k \sigma_i \mathbf{u}_{i}\mathbf{v}_{i}^T = \mathbf{U}_k \mathbf{S}_k \mathbf{V}_k^T.
```

where $\mathbf{U}_k$ and $\mathbf{V}_k$ are the first $k$ columns of $\mathbf{U}$ and $\mathbf{V}$, respectively, and $\mathbf{S}_k$ is the upper-left $k\times k$ submatrix of $\mathbf{S}$.

The rank of a sum of matrices is always less than or equal to the sum of the ranks, so $\mathbf{A}_k$ is a rank-$k$ approximation to $\mathbf{A}$. It turns out that $\mathbf{A}_k$ is the *best* rank-$k$ approximation of $\mathbf{A}$, as measured in the matrix 2-norm.

(thm-best-rank-k)=
::::{proof:theorem}
Suppose $\mathbf{A}$ has rank $r$ and let $\mathbf{A}=\mathbf{U}\mathbf{S}\mathbf{V}^T$ be an SVD. Let $\mathbf{A}_k$ be as in {eq}`svdlowrank` for $1\le k < r$. Then

1. $\| \mathbf{A} - \mathbf{A}_k \|_2 = \sigma_{k+1}, \quad k=1,\ldots,r-1$.
2. If the rank of $\mathbf{B}$ is $k$ or less, then $\| \mathbf{A}-\mathbf{B} \|_2\ge \sigma_{k+1}$.
::::

::::{proof:proof}
(part 1 only) Note that {eq}`svdlowrank` is identical to {eq}`svdsum` with $\sigma_{k+1},\ldots,\sigma_r$ all set to zero. This implies that
  
```{math}
\mathbf{A} - \mathbf{A}_k = \mathbf{U}(\mathbf{S}-\hat{\mathbf{S}})\mathbf{V}^T,
```

where $\hat{\mathbf{S}}$ has those same values of $\sigma_i$ replaced by zero. But that makes the above an SVD of $\mathbf{A} - \mathbf{A}_k$, with singular values $0,\ldots,0,\sigma_{k+1},\ldots,\sigma_r$, the largest of which is $\sigma_{k+1}$. That proves the first claim.
::::

If the singular values of $\mathbf{A}$ decrease sufficiently rapidly, then $\mathbf{A}_{k}$ may capture the most significant behavior of the matrix for a reasonably small value of $k$.


```{index} matrix; as image
```

:::{sidebar} Demo
:class: demo
{doc}`demos/dimreduce-hello`
:::

## Capturing major trends

The use of dimension reduction offered by low-rank SVD approximation goes well beyond simply reducing computation time. By isolating the most important contributions to the matrix,dimension reduction can uncover deep connections and trends that are otherwise obscured by lower-order effects and noise.

One useful way to quantify the decay in the singular values is to compute

```{math}
:label: sing-val-decay
\tau_k = \frac{\sum_{i=1}^k \sigma_i^2}{\sum_{i=1}^r \sigma_i^2}, \quad k=1,\ldots,r.
```

Clearly $0\le \tau_k \le 1$ and $\tau_k$ is non-decreasing as a function of $k$. We can think of $\tau_k$ as the fraction of "energy" contained in the singular values up to and including the $k$th.[^expvar] 

[^expvar]: In statistics this quantity may be interpreted as the fraction of explained variance.

:::{sidebar} Demo
:class: demo
{doc}`demos/dimreduce-voting`
:::

Not all data sets can be reduced effectively to a small number of dimensions, but as {ref}`example-voting` shows, in some cases reduction reveals information that may correspond to real-world understanding.

<!-- \begin{exercises}
	\input{matrixanaly/exercises/DimReduce}
\end{exercises} -->



<!-- \subsection*{Where to learn more}

Details on the computation of the eigenvalue and singular value decompositions are presented at length in~\cite{StewartVol2} and more briefly in Chapters~7 and~8 of~\cite{GolubVan96}. A classic reference on the particulars of the symmetric case is~\cite{Parlett1980}, while~\cite{TrefEmb05} focuses on the non-normal case. Dimension reduction via the SVD often goes by the name *principal component analysis*, which is the subject of~\cite{Jolliffe2002}.
 -->
