---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.12
numbering:
  headings: false
---
# Chapter 7

## Examples

```{code-cell}
:tags: [remove-cell]
include("FNC_init.jl")
```

### 7.1 @section-matrixanaly-insight

(demo-insight-graph-julia)=
``````{dropdown} @demo-insight-graph
:open:
Here we create an adjacency matrix for a graph on four nodes.

```{code-cell}
A = [0 1 0 0; 1 0 0 0; 1 1 0 1; 0 1 1 0]
```

```{index} ! Julia; graphplot
```

The `graphplot` function makes a visual representation of this graph.

```{code-cell}
using Plots, GraphRecipes
graphplot(A, names=1:4, markersize=0.2, arrow=6)
```

Since this adjacency matrix is not symmetric, the edges are all directed, as indicated by the arrows. Here are the counts of all walks of length 3 in the graph:

```{code-cell}
A^3
```

If the adjacency matrix is symmetric, the result is an undirected graph: all edges connect in both directions.

```{code-cell}
A = [0 1 1 0; 1 0 0 1; 1 0 0 0; 0 1 0 0]
graphplot(A, names=1:4, markersize=0.2)
```
``````

(demo-insight-image-julia)=
``````{dropdown} @demo-insight-image
:open:
```{index} ! Julia; Images
```

The `Images` package has many functions for image manipulation, and `TestImages` has some standard images to play with.

```{code-cell}
using Images, TestImages
img = testimage("mandrill")
```

The variable `img` is a matrix.

```{code-cell}
size(img)
```

However, its entries are colors, not numbers.

```{code-cell}
img[100, 10]
```

```{index} ! Julia; eltype
```

You can use `eltype` to find out the type of the elements of any array.

```{code-cell}
eltype(img)
```

It's possible to extract matrices of red, green, and blue intensities, scaled from 0 to 1.

```{code-cell}
R, G, B = red.(img), green.(img), blue.(img);
@show minB, maxB = extrema(B);
```

Or we can convert the pixels to gray, each pixel again scaled from 0 to 1.

```{code-cell}
Gray.(img)
```

In order to do our usual operations, we need to tell Julia that we want to interpret the elements of the image matrix as floating-point values.

```{code-cell}
A = Float64.(Gray.(img))
A[1:4, 1:5]
```

We can use `Gray` to reinterpret a matrix of floating-point values as grayscale pixels.

```{code-cell}
Gray.(reverse(A, dims=1))
```
``````

### 7.2 @section-matrixanaly-evd

(demo-evd-eigen-julia)=
``````{dropdown} @demo-evd-eigen
:open:

```{index} ! Julia; eigvals
```

The `eigvals` function returns a vector of the eigenvalues of a matrix.

```{code-cell}
A = π * ones(2, 2)
```

```{code-cell}
λ = eigvals(A)
```

```{index} ! Julia; eigen
```

If you want the eigenvectors as well, use `eigen`.

```{code-cell}
λ, V = eigen(A)
```

```{code-cell}
norm(A * V[:, 2] - λ[2] * V[:, 2])
```

```{index} ! Julia; sortby
```

Both functions allow you to sort the eigenvalues by specified criteria.

```{code-cell}
A = diagm(-2.3:1.7)
@show eigvals(A, sortby=real);
@show eigvals(A, sortby=abs);
```

If the matrix is not diagonalizable, no message is given, but `V` will be singular. The robust way to detect that circumstance is via $\kappa(\mathbf{V})$.

```{index} condition number; of a matrix
```

```{code-cell}
A = [-1 1; 0 -1]
λ, V = eigen(A)
```

```{code-cell}
cond(V)
```

Even in the nondiagonalizable case, $\mathbf{A}\mathbf{V} = \mathbf{V}\mathbf{D}$ holds.

```{code-cell}
opnorm(A * V - V * diagm(λ))
```
``````

(demo-evd-bauerfike-julia)=
``````{dropdown} @demo-evd-bauerfike
:open:

```{index} Julia; adjoint, Julia; \'
```

We first define a hermitian matrix. Note that the `'` operation is the adjoint and includes complex conjugation.

```{code-cell}
n = 7
A = randn(n, n) + 1im * randn(n, n)
A = (A + A') / 2
```

```{index} Julia; cond
```

We confirm that the matrix $\mathbf{A}$ is normal by checking that $\kappa(\mathbf{V}) = 1$ (to within roundoff).

```{code-cell}
λ, V = eigen(A)
@show cond(V);
```

Now we perturb $\mathbf{A}$ and measure the effect on the eigenvalues. The Bauer–Fike theorem uses absolute differences, not relative ones.
```{tip}
:class: dropdown
Since the ordering of eigenvalues can change, we look at all pairwise differences and take the minima.
```

```{code-cell}
ΔA = 1e-8 * normalize(randn(n, n) + 1im * randn(n, n))
λ̃ = eigvals(A + ΔA)
dist = minimum([abs(x - y) for x in λ̃, y in λ], dims=2)
```

As promised, the perturbations in the eigenvalues do not exceed the normwise perturbation to the original matrix.

Now we see what happens for a triangular matrix.

```{code-cell}
n = 20
x = 1:n
A = triu(x * ones(n)')
A[1:5, 1:5]
```

This matrix is not especially close to normal.

```{code-cell}
λ, V = eigen(A)
@show cond(V);
```

As a result, the eigenvalues can change by a good deal more.

```{code-cell}
ΔA = 1e-8 * normalize(randn(n, n) + 1im * randn(n, n))
λ̃ = eigvals(A + ΔA)
dist = minimum([abs(x - y) for x in λ̃, y in λ], dims=2)
BF_bound = cond(V) * norm(ΔA)
@show maximum(dist), BF_bound;
```

If we plot the eigenvalues of many perturbations, we get a cloud of points that roughly represents all the possible eigenvalues when representing this matrix with single-precision accuracy.

```{code-cell}
using Plots
plt = scatter(λ, zeros(n), aspect_ratio=1)
for _ in 1:200
    ΔA = eps(Float32) * normalize(randn(n, n) + 1im * randn(n, n))
    λ̃ = eigvals(A + ΔA)
    scatter!(real(λ̃), imag(λ̃), m=1, color=:black)
end
plt
```

The plot shows that some eigenvalues are much more affected than others. This situation is not unusual, but it is not explained by the Bauer–Fike theorem.
``````

(demo-evd-francisqr-julia)=
``````{dropdown} @demo-evd-francisqr
:open:
Let's start with a known set of eigenvalues and an orthogonal eigenvector basis.

```{code-cell}
D = diagm([-6, -1, 2, 4, 5])
V, R = qr(randn(5, 5))    # V is unitary
A = V * D * V'
```

```{code-cell}
eigvals(A)
```

Now we will take the QR factorization and just reverse the factors.

```{code-cell}
Q, R = qr(A)
A = R * Q;
```

It turns out that this is a similarity transformation, so the eigenvalues are unchanged.

```{code-cell}
eigvals(A)
```

What's remarkable, and not elementary, is that if we repeat this transformation many times, the resulting matrix converges to $\mathbf{D}$.

```{code-cell}
for k in 1:40
    Q, R = qr(A)
    A = R * Q
end
A
```
``````

### 7.3 @section-matrixanaly-svd

(demo-svd-props-julia)=
``````{dropdown} @demo-svd-props
:open:
We verify some of the fundamental SVD properties using standard Julia functions from `LinearAlgebra`.

```{code-cell}
A = [i^j for i = 1:5, j = 0:3]
```

```{index} ! Julia; svdvals
```


To get only the singular values, use `svdvals`.

```{code-cell}
σ = svdvals(A)
```

Here is verification of the connections between the singular values, norm, and condition number.

```{code-cell}
@show opnorm(A, 2);
@show σ[1];
```

```{code-cell}
@show cond(A, 2);
@show σ[1] / σ[end];
```

```{index} ! Julia; svd
```

To get singular vectors as well, use `svd`. The thin form of the factorization is the default.

```{code-cell}
U, σ, V = svd(A);
@show size(U);
@show size(V);
```

We verify the orthogonality of the singular vectors as follows:

```{code-cell}
@show opnorm(U' * U - I);
@show opnorm(V' * V - I);
```
``````

### 7.4 @section-matrixanaly-symm-eig

(demo-symm-eig-rayleigh-julia)=
``````{dropdown} @demo-symm-eig-rayleigh
:open:
We will use a symmetric matrix with a known EVD and eigenvalues equal to the integers from 1 to 20.

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

If the input to the Rayleigh quotient is within a small $\delta$ of an eigenvector, its output is within $O(\delta^2)$ of the corresponding eigenvalue. In this experiment, we observe that each additional digit of accuracy in an approximate eigenvector gives two more digits to the eigenvalue estimate coming from the Rayleigh quotient.

```{code-cell}
δ = @. 1 ./ 10^(1:5)
eval_diff = zeros(size(δ))
for (k, delta) in enumerate(δ)
    e = randn(n)
    e = delta * e / norm(e)
    x = V[:, 7] + e
    eval_diff[k] = R(x) - 7
end
pretty_table([δ δ.^2 eval_diff];
    column_labels = ["δ (perturbation)", "δ²", "R(x) - λ"], backend=:html)
```
``````

### 7.5 @section-matrixanaly-dimreduce

(demo-dimreduce-hello-julia)=
``````{dropdown} @demo-dimreduce-hello
:open:
We make an image from some text, then reload it as a matrix.

```{code-cell}
using Plots, Images
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
scatter(σ;
    xaxis=(L"i"),  yaxis=(:log10, L"\sigma_i"),
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

(demo-dimreduce-voting-julia)=
``````{dropdown} @demo-dimreduce-voting
:open:
The matrix in [this data file](https://raw.github.com/fncbook/fnc/master/julia/voting.jld2) describes the votes on bills in the 111th session of the United States Senate. (The data set was obtained from [voteview.com](https://voteview.com).) Each row is one senator, and each column is a vote item.

```{code-cell}
using JLD2
@load "voting.jld2" A;
```

If we visualize the votes (yellow is "yea," blue is "nay"), we can see great similarity between many rows, reflecting party unity.

```{code-cell}
heatmap(A;
    color=:viridis,  xlabel="bill",  ylabel="senator",
    title="Votes in 111th U.S. Senate")
```

We use {eq}`sing-val-decay` to quantify the decay rate of the values.

```{code-cell}
U, σ, V = svd(A)
τ = cumsum(σ .^ 2) / sum(σ .^ 2)
scatter(τ[1:16];
    xaxis=("k"),  yaxis=(L"\tau_k"),
    title="Fraction of singular value energy")
```

The first and second singular triples contain about 58% and 17%, respectively, of the energy of the matrix. All others have far less effect, suggesting that the information is primarily two-dimensional. The first left and right singular vectors also contain interesting structure.

```{code-cell}
scatter(U[:, 1], label="", layout=(1, 2),
    xlabel="senator",  title="left singular vector")
scatter!(V[:, 1], label="", subplot=2,
    xlabel="bill",  title="right singular vector")
```

Both vectors have values greatly clustered near $\pm C$ for a constant $C$. These can be roughly interpreted as how partisan a particular senator or bill was, and for which political party. Projecting the senators' vectors into the first two $\mathbf{V}$-coordinates gives a particularly nice way to reduce them to two dimensions. Political scientists label these dimensions *partisanship* and *bipartisanship*. Here we color them by actual party affiliation (also given in the data file): red for Republican, blue for Democrat, and yellow for independent.

```{code-cell}
x1 = A * V[:, 1];
x2 = A * V[:, 2];

@load "voting.jld2" Rep Dem Ind
Rep = vec(Rep);
Dem = vec(Dem);
Ind = vec(Ind);
scatter(x1[Dem], x2[Dem];
    color=:blue,  label="D",
    xaxis=("partisanship"),  yaxis=("bipartisanship"), 
    title="111th US Senate by voting record")
scatter!(x1[Rep], x2[Rep], color=:red, label="R")
scatter!(x1[Ind], x2[Ind], color=:yellow, label="I")
```
``````