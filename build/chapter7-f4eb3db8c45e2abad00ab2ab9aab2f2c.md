---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
# Chapter 7

## Examples

```{code-cell} ipython3
from numpy import *
from numpy.linalg import norm
from matplotlib.pyplot import *
from prettytable import PrettyTable
import sys
sys.path.append('pkg/')
import FNC
import importlib
importlib.reload(FNC)
```

```{code-cell} ipython3
:tags: [remove-cell]
# This (optional) block is for improving the display of plots.
rcParams["figure.figsize"] = [7, 4]
rcParams["lines.linewidth"] = 2
rcParams["lines.markersize"] = 4
rcParams['animation.html'] = "jshtml"  # or try "html5"
``` 

### Section 7.1

(demo-insight-graph-python)=
``````{dropdown} Adjacency matrix
Here we create an adjacency matrix for a graph on four nodes.

```{code-cell}
A = array([[0, 1, 0, 0], [1, 0, 0, 0], [1, 1, 0, 1], [0, 1, 1, 0]])
print(A)
```

The `networkx` package has many functions for working with graphs. Here, we instruct it to create a directed graph from the adjacency matrix, then make a drawing of it.

```{index} ! Python; networkx
```

```{code-cell}
import networkx as nx
G = nx.from_numpy_array(A, create_using=nx.DiGraph)
nx.draw(G, with_labels=True, node_color="yellow")
```

Here are the counts of all walks of length 3 in the graph:

```{code-cell}
print(A**3)
```

If the adjacency matrix is symmetric, the result is an undirected graph: all edges connect in both directions.

```{code-cell}
A = array([[0, 1, 1, 0], [1, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
G = nx.from_numpy_array(A, create_using=nx.Graph)
nx.draw(G, with_labels=True, node_color="yellow")
```
``````

(demo-insight-image-python)=
``````{dropdown} Images as matrices
```{index} ! Julia; Images
```

We will use a test image from the well-known `scikit-image` package.

```{code-cell}
from skimage import data as testimages
img = getattr(testimages, "coffee")()
imshow(img)
```

The variable `img` is a matrix.

```{code-cell}
size(img)
```

However, its entries are colors, not numbers.

```{code-cell}
print(f"image has shape {img.shape}")
print(f"first pixel has value {img[0, 0]}")
```

The three values at each pixel are for intensities of red, green, and blue. We can convert each of those layers into an ordinary matrix of values between 0 and 255, which is maximum intensity.

```{code-cell}
R = img[:, :, 0]
print("upper left corner of the red plane is:")
print(R[:5, :5])
print(f"red channel values range from {R.min()} to {R.max()}")
```

It may also be convenient to convert the image to grayscale, which has just one layer of values from zero (black) to one (white).

```{code-cell}
from skimage.color import rgb2gray
A = rgb2gray(img)
A[:5, :5]
print("upper left corner of grayscale:")
print(A[:5, :5])
print(f"gray values range from {A.min()} to {A.max()}")
```

```{code-cell}
imshow(A, cmap='gray')
axis('off');
```

Some changes we make to the grayscale matrix are easy to interpret visually.

```{code-cell}
imshow(flipud(A), cmap='gray')
axis('off');
```
``````

### Section 7.2

(demo-evd-eigen-python)=
``````{dropdown} Eigenvalues and eigenvectors

```{index} ! Python; eig
```

The `eig` function from `scipy.linalg` will return a vector of eigenvalues and a matrix of associated eigenvectors.

```{code-cell}
from numpy.linalg import eig
A = pi * ones([2, 2])
d, V = eig(A)
print("eigenvalues:", d)
```

We can check the fact that this is an EVD (although in practice we never invert a matrix).

```{code-cell}
from numpy.linalg import inv
D = diag(d)
print(f"should be near zero: {norm(A - V @ D @ inv(V), 2):.2e}")
```

If the matrix is not diagonalizable, no message is given, but `V` will be singular. The robust way to detect that circumstance is via $\kappa(\mathbf{V})$.

```{index} condition number; of a matrix
```

```{index} Python; cond
```

```{code-cell}
from numpy.linalg import cond
A = array([[1, 1], [0, 1]])
d, V = eig(A)
print(f"cond(V) is {cond(V):.2e}")
```

But even in the nondiagonalizable case, $\mathbf{A}\mathbf{V} = \mathbf{V}\mathbf{D}$ holds up to roundoff error.

```{code-cell}
print(f"should be near zero: {norm(A @ V - V @ diag(d), 2):.2e}")
```
```
``````

(demo-evd-bauerfike-python)=
``````{dropdown} Eigenvalue conditioning

We first define a hermitian matrix. Note that we add the *conjugate* transpose of a matrix to itself.

```{code-cell}
n = 7
A = random.randn(n, n) + 1j * random.randn(n, n)
A = (A + conj(A.T)) / 2
```

```{index} Python; cond
```

We confirm that the matrix $\mathbf{A}$ is normal by checking that $\kappa(\mathbf{V}) = 1$ (to within roundoff).

```{code-cell}
from numpy.linalg import eig
d, V = eig(A)
print(f"eigenvector matrix has condition number {cond(V):.5f}")
```

Now we perturb $\mathbf{A}$ and measure the effect on the eigenvalues. Note that the Bauer–Fike theorem uses absolute differences, not relative ones. Since the ordering of eigenvalues can change, we look at all pairwise differences and take the minima.

```{code-cell}
E = random.randn(n, n) + 1j * random.randn(n, n)
E = 1e-8 * E / norm(E, 2)
dd, _ = eig(A + E)
dist = array([min([abs(x - y) for x in dd]) for y in d])
print(dist)
```
As promised, the perturbations in the eigenvalues do not exceed the normwise perturbation to the original matrix.

Now we see what happens for a triangular matrix.

```{code-cell}
n = 20
x = arange(n) + 1
A = triu(outer(x, ones(n)))
print(A[:5, :5])
```

This matrix is not at all close to normal.

```{code-cell}
d, V = eig(A)
print(f"eigenvector matrix has condition number {cond(V):.2e}")
```

As a result, the eigenvalues can change by a good deal more.

```{code-cell}
E = random.randn(n, n) + 1j * random.randn(n, n)
E = 1e-8 * E / norm(E, 2)
dd, _ = eig(A + E)
dist = array([min([abs(x - y) for x in dd]) for y in d])
print(f"Maximum eigenvalue change is {max(dist):.2e}")
print(f"The Bauer-Fike upper bound is {cond(V) * norm(E, 2):.2e}")
```

If we plot the eigenvalues of many perturbations, we get a cloud of points that roughly represents all the possible eigenvalues when representing this matrix with single-precision accuracy.

```{code-cell}
clf
scatter(d, zeros(n), 18)
axis("equal") 
for _ in range(100):
    E = random.randn(n, n) + 1j * random.randn(n, n)
    E = 2**(-26) * E / norm(E, 2)
    dd, _ = eig(A + E)
    scatter(real(dd), imag(dd), 2, 'k')
```

The plot shows that some eigenvalues are much more affected than others. This situation is not unusual, but it is not explained by the Bauer–Fike theorem.
``````

(demo-evd-francisqr-python)=
``````{dropdown} Francis QR iteration
Let's start with a known set of eigenvalues and an orthogonal eigenvector basis.

```{code-cell}
from numpy.linalg import qr
D = diag([-6, -1, 2, 4, 5])
V, R = qr(random.randn(5, 5))
A = V @ D @ V.T    # note that V.T = inv(V) here
```

```{code-cell}
print(sort(eig(A)[0]))
```

Now we will take the QR factorization and just reverse the factors.

```{code-cell}
Q, R = qr(A)
A = R @ Q;
```

It turns out that this is a similarity transformation, so the eigenvalues are unchanged.

```{code-cell}
print(sort(eig(A)[0]))
```

What's remarkable, and not elementary, is that if we repeat this transformation many times, the resulting matrix converges to $\mathbf{D}$.

```{code-cell}
for k in range(40):
    Q, R = qr(A)
    A = R @ Q
set_printoptions(precision=4)
print(A)
```
``````

### Section 7.3

(demo-svd-props-python)=
``````{dropdown} SVD properties
We verify some of the fundamental SVD properties using standard Julia functions from `LinearAlgebra`.

```{code-cell}
A = array([[(i + 1.0) ** j for j in range(4)] for i in range(5)])
set_printoptions(precision=4)
print(A)
```

```{index} ! Python; svd
```

The factorization is obtained using `svd` from `numpy.linalg`.

```{code-cell}
from numpy.linalg import svd
U, sigma, Vh = svd(A)
print("singular values:")
print(sigma)
```

By default, the full factorization type is returned. This can be a memory hog if one of the dimensions of $\mathbf{A}$ is very large.

```{code-cell}
print("size of U:", U.shape)
print("size of V:", Vh.T.shape)
```

Both $\mathbf{U}$ and $\mathbf{V}$ are orthogonal (in the complex case, unitary). Note that it's $\mathbf{V}^*$ that is returned, not $\mathbf{V}$.

```{code-cell}
print(f"should be near zero: {norm(U.T @ U - eye(5), 2):.2e}")
print(f"should be near zero: {norm(Vh @ Vh.T - eye(4), 2):.2e}")
```

Next we test that we have the factorization promised by the SVD, using `diagsvd` to construct a rectangular diagonal matrix.

```{code-cell}
from scipy.linalg import diagsvd
S = diagsvd(sigma, 5, 4)
print(f"should be near zero: {norm(A - U @ S @ Vh, 2):.2e}")
```

Here is verification of the connections between the singular values, norm, and condition number.

```{code-cell}
from numpy.linalg import cond
print("largest singular value:", sigma[0])
print("2-norm of the matrix:  ", norm(A, 2))
print("singular value ratio:", sigma[0] / sigma[-1])
print("2-norm condition no.:", cond(A, 2))
```

For matrices that are much taller than they are wide, the thin SVD form is more memory-efficient.

```{code-cell}
A = random.randn(1000, 10)
U, sigma, Vh = svd(A, full_matrices=False)
print("size of U:", U.shape)
print("size of V:", Vh.shape)
```

``````

### Section 7.4

(demo-symm-eig-rayleigh-python)=
``````{dropdown} Rayleigh quotient
We will use a symmetric matrix with a known EVD.

```{code-cell}
n = 20;
λ = 1:n
D = diagm(λ)
V, _ = qr(randn(n, n))   # get a random orthogonal V
A = V * D * V';
```

The Rayleigh quotient is a scalar-valued function of a vector.

```{code-cell}
R = x -> (x' * A * x) / (x' * x);
```

The Rayleigh quotient evaluated at an eigenvector gives the corresponding eigenvalue.

```{code-cell}
R(V[:, 7])
```

If the input to he Rayleigh quotient is within a small $\delta$ of an eigenvector, its output is within $O(\delta^2)$ of the corresponding eigenvalue. In this experiment, we observe that each additional digit of accuracy in an approximate eigenvector gives two more digits to the eigenvalue estimate coming from the Rayleigh quotient.

```{code-cell}
δ = @. 1 ./ 10^(1:5)
eval_diff = zeros(size(δ))
for (k, delta) in enumerate(δ)
    e = randn(n)
    e = delta * e / norm(e)
    x = V[:, 7] + e
    eval_diff[k] = R(x) - 7
end
labels = ["perturbation δ", "δ²", "R(x) - λ"]
pretty_table([δ δ .^ 2 eval_diff], header=labels)
```
``````

### Section 7.5

(demo-dimreduce-hello-python)=
``````{dropdown} Image compression
We make an image from some text, then reload it as a matrix.

```{code-cell}
plot(annotations=(0.5, 0.5, text("Hello world", 44, :center, :center)),
    grid=:none, frame=:none, size=(400, 150))
savefig("hello.png")
img = load("hello.png")
A = @. Float64(Gray(img))
Gray.(A)
```

Next we show that the singular values decrease until they reach zero (more precisely, until they are about $\epsilon_\text{mach}$ times the norm of the matrix) at around $k=45$.

```{code-cell}
U, σ, V = svd(A)
scatter(σ, xaxis=(L"i"), yaxis=(:log10, L"\sigma_i"),
    title="Singular values")
```

The rapid decrease suggests that we can get fairly good low-rank approximations.

```{code-cell}
plt = plot(layout=(2, 2), frame=:none, aspect_ratio=1, titlefontsize=10)
for i in 1:4
    k = 3i
    Ak = U[:, 1:k] * diagm(σ[1:k]) * V[:, 1:k]'
    plot!(Gray.(Ak), subplot=i, title="rank = $k")
end
plt
```

Consider how little data is needed to reconstruct these images. For rank-9, for instance, we have 9 left and right singular vectors plus 9 singular values, for a compression ratio of better than 12:1.

```{code-cell}
m, n = size(A)
compression = m * n / (9 * (m + n + 1))
```
``````

(demo-dimreduce-voting-python)=
``````{dropdown} Dimension reduction in voting records
This matrix describes the votes on bills in the 111th session of the United States Senate. (The data set was obtained from [https://voteview.com].) Each row is one senator, and each column is a vote item.

```{code-cell}
@load "voting.jld2" A;
```

If we visualize the votes (yellow is "yea," blue is "nay"), we can see great similarity between many rows, reflecting party unity.

```{code-cell}
heatmap(A, color=:viridis,
    title="Votes in 111th U.S. Senate", xlabel="bill", ylabel="senator")
```

We use {eq}`sing-val-decay` to quantify the decay rate of the values.

```{code-cell}
U, σ, V = svd(A)
τ = cumsum(σ .^ 2) / sum(σ .^ 2)
scatter(τ[1:16], xaxis=("k"), yaxis=(L"\tau_k"),
    title="Fraction of singular value energy")
```

The first and second singular triples contain about 58% and 17%, respectively, of the energy of the matrix. All others have far less effect, suggesting that the information is primarily two-dimensional. The first left and right singular vectors also contain interesting structure.

```{code-cell}
scatter(U[:, 1], label="", layout=(1, 2),
    xlabel="senator", title="left singular vector")
scatter!(V[:, 1], label="", subplot=2,
    xlabel="bill", title="right singular vector")
```

Both vectors have values greatly clustered near $\pm C$ for a constant $C$. These can be roughly interpreted as how partisan a particular senator or bill was, and for which political party. Projecting the senators' vectors into the first two $\mathbf{V}$-coordinates gives a particularly nice way to reduce them to two dimensions. Political scientists label these dimensions *partisanship* and *bipartisanship*. Here we color them by actual party affiliation (also given in the data file): red for Republican, blue for Democrat, and yellow for independent.

```{code-cell}
x1 = A * V[:, 1];
x2 = A * V[:, 2];

@load "voting.jld2" Rep Dem Ind
Rep = vec(Rep);
Dem = vec(Dem);
Ind = vec(Ind);
scatter(x1[Dem], x2[Dem], color=:blue, label="D",
    xaxis=("partisanship"), yaxis=("bipartisanship"), title="111th US Senate by voting record")
scatter!(x1[Rep], x2[Rep], color=:red, label="R")
scatter!(x1[Ind], x2[Ind], color=:yellow, label="I")
```
``````