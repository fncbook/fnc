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
exec(open("FNC_init.py").read())
```

### 7.1 @section-matrixanaly-insight

(demo-insight-graph-python)=
``````{dropdown} @demo-insight-graph
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
``````{dropdown} @demo-insight-image
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

### 7.2 @section-matrixanaly-evd

(demo-evd-eigen-python)=
``````{dropdown} @demo-evd-eigen

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
``````{dropdown} @demo-evd-bauerfike

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
    E = finfo(np.float32).eps * E / norm(E, 2)
    dd, _ = eig(A + E)
    scatter(real(dd), imag(dd), 2, 'k')
```

The plot shows that some eigenvalues are much more affected than others. This situation is not unusual, but it is not explained by the Bauer–Fike theorem.
``````

(demo-evd-francisqr-python)=
``````{dropdown} @demo-evd-francisqr
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

### 7.3 @section-matrixanaly-svd

(demo-svd-props-python)=
``````{dropdown} @demo-svd-props
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

For matrices that are much taller than they are wide, the thin SVD form is more memory-efficient, because $\mathbf{U}$ takes the same shape.

```{code-cell}
A = random.randn(1000, 10)
U, sigma, Vh = svd(A, full_matrices=False)
print("size of U:", U.shape)
print("size of V:", Vh.shape)
```
``````

### 7.4 @section-matrixanaly-symm-eig

(demo-symm-eig-rayleigh-python)=
``````{dropdown} @demo-symm-eig-rayleigh
We will use a symmetric matrix with a known EVD and eigenvalues equal to the integers from 1 to 20.

```{code-cell}
from numpy.linalg import qr
n = 20
d = arange(n) + 1
D = diag(d)
V, _ = qr(random.randn(n, n))    # get a random orthogonal V
A = V @ D @ V.T
```

The Rayleigh quotient is a scalar-valued function of a vector.

```{code-cell}
R = lambda x: dot(x, A @ x) / dot(x, x)
```

The Rayleigh quotient evaluated at an eigenvector gives the corresponding eigenvalue.

```{code-cell}
print(R(V[:, 6]))
```

If the input to he Rayleigh quotient is within a small $\delta$ of an eigenvector, its output is within $O(\delta^2)$ of the corresponding eigenvalue. In this experiment, we observe that each additional digit of accuracy in an approximate eigenvector gives two more digits to the eigenvalue estimate coming from the Rayleigh quotient.

```{code-cell}
results = PrettyTable(["perturbation size", "R.Q. - λ"])
for delta in 1 / 10 ** arange(1, 6):
    e = random.randn(n)
    e = delta * e / norm(e)
    x = V[:, 5] + e
    quotient = R(x)
    results.add_row([delta, quotient - d[5]])

print(results)
```
``````

### 7.5 @section-matrixanaly-dimreduce

(demo-dimreduce-hello-python)=
``````{dropdown} @demo-dimreduce-hello
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
``````

(demo-dimreduce-voting-python)=
``````{dropdown} @demo-dimreduce-voting
This matrix describes the votes on bills in the 111th session of the United States Senate. (The data set was obtained from [https://voteview.com].) Each row is one senator, and each column is a vote item.

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
``````
