---
numbering:
  enumerator: 7.5.%s
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
text(
    0.5,
    0.5,
    "Hello world",
    dict(fontsize=44),
    horizontalalignment="center",
    verticalalignment="center",
)
axis("off")
savefig("hello.png")
```

```{code-cell}
from skimage.io import imread
from skimage.color import rgb2gray
img = imread("hello.png")[:, :, :3]
A = rgb2gray(img)
print(f"image of size {A.shape}")
```


Next we show that the singular values decrease until they reach zero (more precisely, until they are about $\epsilon_\text{mach}$ times the norm of the matrix) at around index $38$.

```{code-cell}
from numpy.linalg import svd
U, sigma, Vt = svd(A)
semilogy(sigma, "o")
title("Singular values")
xlabel("$i$"), ylabel("$\\sigma_i$");
```

```{code-cell}
significant = sigma / sigma[0] > 10 * 2**-52
print(f"last significant singular value at index {max(where(significant)[0])}")
```

The rapid decrease suggests that we can get fairly good low-rank approximations.

```{code-cell}
for k in range(4):
    r = 2 + 2 * k
    Ak = U[:, :r] @ diag(sigma[:r]) @ Vt[:r, :]
    subplot(2, 2, k + 1)
    imshow(Ak, cmap="gray", clim=(0.0, 1.0))
    title(f"rank = {r}")
    xticks([]), yticks([])
```

Consider how little data is needed to reconstruct these images. For rank-8, for instance, we have 8 left and right singular vectors plus 8 singular values.

```{code-cell}
m, n = A.shape
full_size = m * n
compressed_size = 8 * (m + n + 1)
print(f"compression ratio: {full_size / compressed_size:.1f}")
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

The matrix in [this data file](https://raw.github.com/fncbook/fnc/master/python/voting.mat) describes the votes on bills in the 111th session of the United States Senate. (The data set was obtained from [voteview.com](https://voteview.com).) Each row is one senator, and each column is a vote item.

```{code-cell}
from scipy.io import loadmat
vars = loadmat("voting.mat")
A = vars["A"]
m, n = A.shape
print("size:", (m, n))
```

If we visualize the votes (yellow is "yea," blue is "nay"), we can see great similarity between many rows, reflecting party unity.

```{code-cell}
imshow(A, cmap="viridis")
xlabel("bill")
ylabel("senator")
title("Votes in 111th U.S. Senate");
```

We use {eq}`sing-val-decay` to quantify the decay rate of the values.

```{code-cell}
U, sigma, Vt = svd(A)
tau = cumsum(sigma**2) / sum(sigma**2)
plot(range(1, 17), tau[:16], "o")
xlabel("$k$")
ylabel("$\tau_k$")
title("Fraction of singular value energy");
```

The first and second singular triples contain about 58% and 17%, respectively, of the energy of the matrix. All others have far less effect, suggesting that the information is primarily two-dimensional. The first left and right singular vectors also contain interesting structure.

```{code-cell}
subplot(1, 2, 1)
plot(U[:, 0], "o")
xlabel("senator"),title("left singular vector")
subplot(1, 2, 2)
plot(Vt[0, :], "o")
xlabel("bill"), title("right singular vector");
```

Both vectors have values greatly clustered near $\pm C$ for a constant $C$. These can be roughly interpreted as how partisan a particular senator or bill was, and for which political party. Projecting the senators' vectors into the first two $\mathbf{V}$-coordinates gives a particularly nice way to reduce them to two dimensions. Political scientists label these dimensions *partisanship* and *bipartisanship*. Here we color them by actual party affiliation (also given in the data file): red for Republican, blue for Democrat, and yellow for independent.

```{code-cell}
x1 = sigma[0] * U[:, 0]
x2 = sigma[1] * U[:, 1]

Rep = vars["Rep"] - 1
Dem = vars["Dem"] - 1
Ind = vars["Ind"] - 1

scatter(x1[Dem], x2[Dem], color="blue", label="D")
scatter(x1[Rep], x2[Rep], color="red", label="R")
scatter(x1[Ind], x2[Ind], color="darkorange", label="I")

xlabel("partisanship"),  ylabel("bipartisanship")
legend(),  title("111th US Senate in 2D");
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
