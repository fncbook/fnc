---
numbering:
  enumerator: 7.5.%s
---
(section-matrixanaly-dimreduce)=
# Dimension reduction

```{index} singular value decomposition
```
```{index} dimension reduction
```

The SVD has another important property that proves very useful in a variety of applications. Let $\mathbf{A}$ be a real $m\times n$ matrix with SVD $\mathbf{A}=\mathbf{U}\mathbf{S}\mathbf{V}^T$ and (momentarily) $m\ge n$. Another way of writing the thin form of the SVD is

::::{math}
:label: svdsum
\begin{split}
  \mathbf{A} = \hat{\mathbf{U}}\hat{\mathbf{S}}\mathbf{V}^T &=
  \begin{bmatrix}
    \rule[-0.3em]{0pt}{1em} \mathbf{u}_1 & \mathbf{u}_2 & \cdots & \mathbf{u}_n
  \end{bmatrix} \:
  \begin{bmatrix}
    \sigma_1 & & \\
    & \ddots & \\
    & & \sigma_n
  \end{bmatrix} \: 
        \begin{bmatrix}
          \mathbf{v}_1^T \\ \vdots \\ \mathbf{v}_n^T
        \end{bmatrix}\ \\
  &=
  \begin{bmatrix}
    \rule[-0.3em]{0pt}{1em} \sigma_1\mathbf{u}_1  & \cdots & \sigma_n\mathbf{u}_n
  \end{bmatrix}\:
  \begin{bmatrix}
    \mathbf{v}_1^T \\ \vdots \\ \mathbf{v}_n^T
  \end{bmatrix} \\
  &= \sigma_1 \mathbf{u}_{1}\mathbf{v}_{1}^T + \cdots + \sigma_r \mathbf{u}_{r}\mathbf{v}_{r}^T = \sum_{i=1}^r \sigma_i \mathbf{u}_{i}\mathbf{v}_{i}^T,
\end{split}
::::

where $r$ is the rank of $\mathbf{A}$. The final formula also holds for the case $m<n$.

Each outer product $\mathbf{u}_{i}\mathbf{v}_{i}^T$ is a rank-1 matrix of unit 2-norm. Thanks to the ordering of singular values, then, Equation {eq}`svdsum` expresses $\mathbf{A}$ as a sum of decreasingly important contributions. This motivates the definition, for $1\le k \le r$,

```{math}
:label: svdlowrank
\mathbf{A}_k = \sum_{i=1}^k \sigma_i \mathbf{u}_{i}\mathbf{v}_{i}^T = \mathbf{U}_k \mathbf{S}_k \mathbf{V}_k^T,
```

where $\mathbf{U}_k$ and $\mathbf{V}_k$ are the first $k$ columns of $\mathbf{U}$ and $\mathbf{V}$, respectively, and $\mathbf{S}_k$ is the upper-left $k\times k$ submatrix of $\mathbf{S}$.

The rank of a sum of matrices is always less than or equal to the sum of the ranks, so $\mathbf{A}_k$ is a rank-$k$ approximation to $\mathbf{A}$. It turns out that $\mathbf{A}_k$ is the *best* rank-$k$ approximation of $\mathbf{A}$, as measured in the matrix 2-norm.

(theorem-best-rank-k)=
::::{prf:theorem}
Suppose $\mathbf{A}$ has rank $r$ and let $\mathbf{A}=\mathbf{U}\mathbf{S}\mathbf{V}^T$ be an SVD. Let $\mathbf{A}_k$ be as in {eq}`svdlowrank` for $1\le k < r$. Then

1. $\| \mathbf{A} - \mathbf{A}_k \|_2 = \sigma_{k+1}, \quad k=1,\ldots,r-1$, and
2. If the rank of $\mathbf{B}$ is $k$ or less, then $\| \mathbf{A}-\mathbf{B} \|_2\ge \sigma_{k+1}$.
::::

::::{prf:proof}
(part 1 only) Note that {eq}`svdlowrank` is identical to {eq}`svdsum` with $\sigma_{k+1},\ldots,\sigma_r$ all set to zero. This implies that
  
```{math}
\mathbf{A} - \mathbf{A}_k = \mathbf{U}(\mathbf{S}-\hat{\mathbf{S}})\mathbf{V}^T,
```

where $\hat{\mathbf{S}}$ has those same values of $\sigma_i$ replaced by zero. But that makes the above an SVD of $\mathbf{A} - \mathbf{A}_k$, with singular values $0,\ldots,0,\sigma_{k+1},\ldots,\sigma_r$, the largest of which is $\sigma_{k+1}$. That proves the first claim.
::::

## Compression

If the singular values of $\mathbf{A}$ decrease sufficiently rapidly, then $\mathbf{A}_{k}$ may capture the most significant behavior of the matrix for a reasonably small value of $k$.


```{index} image (as a matrix)
```

(demo-dimreduce-hello)=
```````{prf:example}
``````{tab-set}
`````{tab-item} Julia
:sync: julia
````{dropdown} Image compression
:open:
```{include} julia/hello.ipynb
```
````
`````

`````{tab-item} MATLAB
:sync: matlab
````{dropdown} Image compression
:open:
```{include} matlab/hello.ipynb
```
````
`````

`````{tab-item} Python
:sync: python
````{dropdown} Image compression
:open:
```{include} python/hello.ipynb
```
````
`````
``````
```````
    

## Capturing major trends

The use of dimension reduction offered by low-rank SVD approximation goes well beyond simply reducing computation time. By isolating the most important contributions to the matrix, dimension reduction can uncover deep connections and trends that are otherwise obscured by weaker effects and noise.

One useful way to quantify the decay in the singular values is to compute

```{math}
:label: sing-val-decay
s_k = \sum_{i=1}^k \sigma_i^2, \quad \tau_k = \frac{s_k}{s_r}, \quad k=1,\ldots,r.
```

Clearly $0\le \tau_k \le 1$ and $\tau_k$ is non-decreasing as a function of $k$. We can think of $\tau_k$ as the fraction of energy (or in statistical terms, variance) contained in the singular values up to and including the $k$th.[^expvar] 

[^expvar]: In statistics this quantity may be interpreted as the fraction of explained variance.

(demo-dimreduce-voting)=
```````{prf:example}
``````{tab-set}
`````{tab-item} Julia
:sync: julia
````{dropdown} Dimension reduction in voting records
:open:
```{include} julia/voting.ipynb
```
````
`````

`````{tab-item} MATLAB
:sync: matlab
````{dropdown} Dimension reduction in voting records
:open:
```{include} matlab/voting.ipynb
```
````
`````

`````{tab-item} Python
:sync: python
````{dropdown} Dimension reduction in voting records
:open:
```{include} python/voting.ipynb
```
````
`````
``````
```````
    

Not all data sets can be reduced effectively to a small number of dimensions, but as {numref}`Demo {number} <demo-dimreduce-voting>` illustrates, in some cases reduction reveals information that corresponds to real-world understanding.

## Exercises

1. ✍  Suppose that $\mathbf{A}$ is an $n\times n$ matrix. Explain why $\sigma_n$ is the distance (in 2-norm) from $\mathbf{A}$ to the set of all singular matrices.

2. ✍ Suppose $\mathbf{A}$ is a $7\times 4$ matrix and the eigenvalues of $\mathbf{A}^*\mathbf{A}$ are 3, 4, 7, and 10. How close is $\mathbf{A}$ in the 2-norm to (a) a rank-3 matrix? (b) a rank-2 matrix? 

3. 
    **(a)** ⌨ Find the rank-1 matrix closest to 
    
    $$
    \mathbf{A}=\displaystyle \begin{bmatrix}
    1 & 5 \\ 5 & 1
    \end{bmatrix},
    $$
    
    as measured in the 2-norm.

    **(b)** ⌨ Repeat part (a) for 
    
    $$
    \mathbf{A}=\displaystyle \begin{bmatrix}
    1 & 5 \\ 0 & 1
    \end{bmatrix}.
    $$

4. ✍ Find the rank-1 matrix closest to 
   
    $$
    \mathbf{A}=\displaystyle \begin{bmatrix}
    1 & b \\ b & 1
    \end{bmatrix},
    $$
  
    as measured in the 2-norm, where $b>0$.

5. ⌨ Following {numref}`Demo %s <demo-dimreduce-hello>` as a guide, load the "mandrill" test image and convert it to a matrix of floating-point pixel grayscale intensities. Using the SVD, display as images the best approximations of rank 5, 10, 15, and 20. 
