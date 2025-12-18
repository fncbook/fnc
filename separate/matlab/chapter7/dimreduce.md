---
numbering:
  enumerator: 7.5.%s
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
---
```{code-cell}
:tags: [remove-cell]
clear all
format short
set(0, 'defaultaxesfontsize', 12)
set(0, 'defaultlinelinewidth', 1.5)
set(0, 'defaultFunctionLinelinewidth', 1.5)
set(0, 'defaultscattermarkerfacecolor', 'flat')
gcf;
set(gcf, 'Position', [0 0 600 350])
addpath FNC-matlab
```

(section-matrixanaly-dimreduce)=

# Dimension reduction

```{index} singular value decomposition
```

```{index} dimension reduction
```

The SVD has another important property that proves very useful in a variety of applications. Let $\mathbf{A}$ be a real $m\times n$ matrix with SVD $\mathbf{A}=\mathbf{U}\mathbf{S}\mathbf{V}^T$ and (momentarily) $m\ge n$. Another way of writing the thin form of the SVD is

````{math}
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
````

where $r$ is the rank of $\mathbf{A}$. The final formula also holds for the case $m<n$.

Each outer product $\mathbf{u}_{i}\mathbf{v}_{i}^T$ is a rank-1 matrix of unit 2-norm. Thanks to the ordering of singular values, then, Equation {eq}`svdsum` expresses $\mathbf{A}$ as a sum of decreasingly important contributions. This motivates the definition, for $1\le k \le r$,

```{math}
:label: svdlowrank
\mathbf{A}_k = \sum_{i=1}^k \sigma_i \mathbf{u}_{i}\mathbf{v}_{i}^T = \mathbf{U}_k \mathbf{S}_k \mathbf{V}_k^T,
```

where $\mathbf{U}_k$ and $\mathbf{V}_k$ are the first $k$ columns of $\mathbf{U}$ and $\mathbf{V}$, respectively, and $\mathbf{S}_k$ is the upper-left $k\times k$ submatrix of $\mathbf{S}$.

The rank of a sum of matrices is always less than or equal to the sum of the ranks, so $\mathbf{A}_k$ is a rank-$k$ approximation to $\mathbf{A}$. It turns out that $\mathbf{A}_k$ is the *best* rank-$k$ approximation of $\mathbf{A}$, as measured in the matrix 2-norm.

::::{prf:theorem}
:label: theorem-best-rank-k
Suppose $\mathbf{A}$ has rank $r$ and let $\mathbf{A}=\mathbf{U}\mathbf{S}\mathbf{V}^T$ be an SVD. Let $\mathbf{A}_k$ be as in {eq}`svdlowrank` for $1\le k < r$. Then

1. $\| \mathbf{A} - \mathbf{A}_k \|_2 = \sigma_{k+1}, \quad k=1,\ldots,r-1$, and
2. If the rank of $\mathbf{B}$ is $k$ or less, then $\| \mathbf{A}-\mathbf{B} \|_2\ge \sigma_{k+1}$.
::::

::::{prf:proof}
:enumerated: false

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

::::{prf:example} Image compression
:label: demo-dimreduce-hello

We make an image from some text, then reload it as a matrix.

```{code-cell}
:tags: hide-input
clf
tobj = text(0, 0,'Hello world','fontsize',44);
ex = get(tobj, 'extent');
axis([ex(1) ex(1) + ex(3) ex(2) ex(2) + ex(4)]), axis off
exportgraphics(gca, 'hello.png', resolution=300)
A = imread('hello.png');
A = double(im2gray(A));
size_A = size(A)
```

Next we show that the singular values decrease until they reach zero (more precisely, until they are about $\epsilon_\text{mach}$ times the norm of the matrix) at around $k=100$.

```{code-cell}
[U, S, V] = svd(A);
sigma = diag(S);
semilogy(sigma, '.')
title('singular values'), axis tight 
xlabel('i'), ylabel('\sigma_i') 
r = find(sigma / sigma(1) > 10*eps, 1, 'last')
```

The rapid decrease suggests that we can get fairly good low-rank approximations.

```{code-cell}
for i = 1:4
    subplot(2, 2, i)
    k = 2*i;
    Ak = U(:, 1:k) * S(1:k, 1:k) * V(:, 1:k)';
    imshow(Ak, [0, 255])
    title(sprintf('rank = %d', k))
end
```

Consider how little data is needed to reconstruct these images. For rank-9, for instance, we have 9 left and right singular vectors plus 9 singular values, for a compression ratio of better than 12:1.

```{code-cell}
[m, n] = size(A);
full_size = m * n;
compressed_size = 8 * (m + n + 1);
fprintf("compression ratio: %.1f", full_size / compressed_size)
```

::::

## Capturing major trends

The use of dimension reduction offered by low-rank SVD approximation goes well beyond simply reducing computation time. By isolating the most important contributions to the matrix, dimension reduction can uncover deep connections and trends that are otherwise obscured by weaker effects and noise.

One useful way to quantify the decay in the singular values is to compute

```{math}
:label: sing-val-decay
s_k = \sum_{i=1}^k \sigma_i^2, \quad \tau_k = \frac{s_k}{s_r}, \quad k=1,\ldots,r.
```

Clearly $0\le \tau_k \le 1$ and $\tau_k$ is non-decreasing as a function of $k$. We can think of $\tau_k$ as the fraction of energy (or in statistical terms, variance) contained in the singular values up to and including the $k$th.[^expvar]

[^expvar]: In statistics this quantity may be interpreted as the fraction of explained variance.

::::{prf:example} Dimension reduction in voting records
:label: demo-dimreduce-voting

The matrix in [this data file](https://raw.github.com/fncbook/fnc/master/matlab/voting.mat) describes the votes on bills in the 111th session of the United States Senate. (The data set was obtained from [voteview.com](https://voteview.com).) Each row is one senator, and each column is a vote item.

```{code-cell}
load voting
```

If we visualize the votes (yellow is "yea," blue is "nay"), we can see great similarity between many rows, reflecting party unity.

```{code-cell}
clf
imagesc(A)
colormap parula
title('Votes in 111th U.S. Senate')
ylabel(('senator'),  xlabel('bill'));
```

We use {eq}`sing-val-decay` to quantify the decay rate of the values.

```{code-cell}
[U, S, V] = svd(A);
sigma = diag(S);
tau = cumsum(sigma.^2) / sum(sigma.^2);
plot(tau(1:16), 'o')
xlabel('k'),  ylabel('\tau_k')
title(('Fraction of singular value energy'));
```

The first and second singular triples contain about 58% and 17%, respectively, of the energy of the matrix. All others have far less effect, suggesting that the information is primarily two-dimensional. The first left and right singular vectors also contain interesting structure.

```{code-cell}
subplot(211), plot(U(:, 1), '.')
xlabel('senator number'), title('left singular vector')
subplot(212), plot(V(:, 1), '.')
xlabel('bill number'), title(('right singular vector'));
```

Both vectors have values greatly clustered near $\pm C$ for a constant $C$. These can be roughly interpreted as how partisan a particular senator or bill was, and for which political party. Projecting the senators' vectors into the first two $\mathbf{V}$-coordinates gives a particularly nice way to reduce them to two dimensions. Political scientists label these dimensions *partisanship* and *bipartisanship*. Here we color them by actual party affiliation (also given in the data file): red for Republican, blue for Democrat, and yellow for independent.

```{code-cell}
clf
x1 = V(:, 1)'*A';   x2 = V(:, 2)'*A'; 
scatter(x1(Dem), x2(Dem), 20, 'b'),  hold on
scatter(x1(Rep), x2(Rep), 20, 'r')
scatter(x1(Ind), x2(Ind), 20, 'm')
xlabel('partisanship'),  ylabel('bipartisanship')
legend('Democrat', 'Republican', 'Independent')
title(('111th US Senate in 2D'));
```

::::

Not all data sets can be reduced effectively to a small number of dimensions, but as @demo-dimreduce-voting illustrates, in some cases reduction reveals information that corresponds to real-world understanding.

## Exercises

``````{exercise}
:label: problem-dimreduce-distance
✍  Suppose that $\mathbf{A}$ is an $n\times n$ matrix. Explain why $\sigma_n$ is the distance (in 2-norm) from $\mathbf{A}$ to the set of all singular matrices.
``````

``````{exercise}
:label: problem-dimreduce-rank
✍ Suppose $\mathbf{A}$ is a $7\times 4$ matrix and the eigenvalues of $\mathbf{A}^*\mathbf{A}$ are 3, 4, 7, and 10. How close is $\mathbf{A}$ in the 2-norm to (a) a rank-3 matrix? (b) a rank-2 matrix? 
``````

``````{exercise}
:label: problem-dimreduce-rank1
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
``````

``````{exercise}
:label: problem-dimreduce-rank1b
✍ Find the rank-1 matrix closest to 

$$
\mathbf{A}=\displaystyle \begin{bmatrix}
1 & b \\ b & 1
\end{bmatrix},
$$

as measured in the 2-norm, where $b>0$. It may help to know that $\mathbf{A}$ has the EVD

```{math}
:numbered: false
\left(
\frac{1}{2}
\begin{bmatrix}
1 & 1 \\ 1 & -1
\end{bmatrix}
\right)
\begin{bmatrix}
1 + b & 0 \\ 0 & 1 - b
\end{bmatrix}
\begin{bmatrix}
1 & 1 \\ 1 & -1
\end{bmatrix}.
```
``````

``````{exercise}
:label: problem-dimreduce-image
⌨ Following @demo-dimreduce-hello as a guide, load the "mandrill" test image and convert it to a matrix of floating-point pixel grayscale intensities. Using the SVD, display as images the best approximations of rank 5, 10, 15, and 20. 
``````
