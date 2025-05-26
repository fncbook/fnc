---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---
# Chapter 7

## Examples

```{code-cell}
:tags: [remove-cell]
cd  /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init;
pwd;
```

### 7.1 @section-matrixanaly-insight

(demo-insight-graph-matlab)=
``````{dropdown} @demo-insight-graph
:open:
Here we create an adjacency matrix for a graph on four nodes.

```{code-cell}
A = [0 1 0 0; 1 0 0 0; 1 1 0 1; 0 1 1 0]
```

```{index} ! MATLAB; graph (network)
```

Since this adjacency matrix is not symmetric, the edges are all directed. We use `digraph` to create a directed graph.

```{code-cell}
G = digraph(A);
plot(G)
```
Here are the counts of all walks of length 3 in the graph:

```{code-cell}
A^3
```

If the adjacency matrix is symmetric, the result is an undirected graph: all edges connect in both directions.

```{code-cell}
A = [0 1 1 0; 1 0 0 1; 1 0 0 0; 0 1 0 0];
plot(graph(A))
```

A "buckyball" is an allotrope of carbon atoms with the same connection structure as a soccer ball.

```{code-cell}
plot(graph(bucky))
```
``````

(demo-insight-image-matlab)=
``````{dropdown} @demo-insight-image
:open:
```{index} ! Julia; Images
```

MATLAB ships with a few test images to play with.

```{code-cell}
A = imread('peppers.png');
color_size = size(A)
```

Use `imshow` to display the image.

```{code-cell}
imshow(A)
```

The image has three layers or channels for red, green, and blue. We can deal with each layer as a matrix, or (as below) convert it to a single matrix indicating shades of gray from black (0) to white (255). Either way, we have to explicitly convert the entries to floating-point values rather than integers. 


```{code-cell}
A = im2gray(A);   % collapse from 3 dimensions to 2
gray_size = size(A)
imshow(A)
```

Before we can do any numerical computation, we need to convert the image to a matrix of floating-point numbers.

```{code-cell}
A = double(A);
```
``````

### 7.2 @section-matrixanaly-evd

(demo-evd-eigen-matlab)=
``````{dropdown} @demo-evd-eigen
:open:

```{index} ! MATLAB; eig
```

The `eig` function with one output argument returns a vector of the eigenvalues of a matrix.

```{code-cell}
A = pi * ones(2, 2);
lambda = eig(A)
```

With two output arguments given, `eig` returns a matrix eigenvectors and a diagonal matrix with the eigenvalues.

```{code-cell}
[V, D] = eig(A)
```

We can check the fact that this is an EVD.

```{code-cell}
norm( A - V*D/V )   % / V is like * inv(V)
```

If the matrix is not diagonalizable, no message is given, but `V` will be singular. The robust way to detect that circumstance is via $\kappa(\mathbf{V})$.

```{index} condition number; of a matrix
```

```{code-cell}
A = [-1 1; 0 -1];
[V, D] = eig(A)
```

```{code-cell}
cond(V)
```

Even in the nondiagonalizable case, $\mathbf{A}\mathbf{V} = \mathbf{V}\mathbf{D}$ holds.

```{code-cell}
norm(A * V - V * D)
```
``````

(demo-evd-bauerfike-matlab)=
``````{dropdown} @demo-evd-bauerfike
:open:

```{index} MATLAB; adjoint, MATLAB; \'
```

We first define a hermitian matrix. Note that the `'` operation is the adjoint and includes complex conjugation.

```{code-cell}
n = 7;
A = randn(n, n) + 1i * randn(n, n);
A = (A + A') / 2;
```

We confirm that the matrix $\mathbf{A}$ is normal by checking that $\kappa(\mathbf{V}) = 1$ (to within roundoff).

```{code-cell}
[V, D] = eig(A);
lambda = diag(D);
cond(V)
```

Now we perturb $\mathbf{A}$ and measure the effect on the eigenvalues. The Bauer–Fike theorem uses absolute differences, not relative ones. Note: since the ordering of eigenvalues can change, we look at all pairwise differences and take the minima.

```{code-cell}
E = randn(n, n) + 1i * randn(n, n);
E = 1e-8 * E / norm(E);
dd = eig(A + E);
dist = [];
for j = 1:n
    dist = [dist; min(abs(dd - lambda(j)))];
end
dist
```

As promised, the perturbations in the eigenvalues do not exceed the normwise perturbation to the original matrix.

Now we see what happens for a triangular matrix.

```{code-cell}
n = 20;
x = (1:n)';
A = triu(x * ones(1, n));
A(1:5, 1:5)
```

This matrix is not at all close to normal.

```{code-cell}
[V, D] = eig(A);
lambda = diag(D);
cond(V)
```

As a result, the eigenvalues can change by a good deal more.

```{code-cell}
E = randn(n, n) + 1i * randn(n, n);
E = 1e-8 * E / norm(E);
dd = eig(A + E);
dist = -Inf;
for j = 1:n
    dist = max(dist, min(abs(dd - lambda(j))));
end
fprintf("max change in eigenvalues: %.2e", dist)
fprintf("Bauer-Fike upper bound: %.2e", cond(V) * norm(E))
```

If we plot the eigenvalues of many perturbations, we get a cloud of points that roughly represents all the possible eigenvalues when representing this matrix with single-precision accuracy.

```{code-cell}
clf
scatter(lambda, 0*lambda)
axis equal; hold on
for k = 1:60
    E = randn(n, n) + 1i * randn(n, n);
    E = eps(single(1)) * E / norm(E);
    dd = eig(A + E);
    plot(real(dd), imag(dd), 'k.', markersize=2)
end
```

The plot shows that some eigenvalues are much more affected than others. This situation is not unusual, but it is not explained by the Bauer–Fike theorem.
``````

(demo-evd-francisqr-matlab)=
``````{dropdown} @demo-evd-francisqr
:open:
Let's start with a known set of eigenvalues and an orthogonal eigenvector basis.

```{code-cell}
D = diag([-6, -1, 2, 4, 5]);
[V, R]= qr(randn(5, 5));    % V is unitary
A = V * D * V';
```

```{code-cell}
sort(eig(A))
```

Now we will take the QR factorization and just reverse the factors.

```{code-cell}
[Q, R] = qr(A);
A = R * Q;
```

It turns out that this is a similarity transformation, so the eigenvalues are unchanged.

```{code-cell}
sort(eig(A))
```

What's remarkable, and not elementary, is that if we repeat this transformation many times, the resulting matrix converges to $\mathbf{D}$.

```{code-cell}
for k = 1:40
    [Q, R] = qr(A);
    A = R * Q;
end
format short e
A
```
``````

### 7.3 @section-matrixanaly-svd

(demo-svd-props-matlab)=
``````{dropdown} @demo-svd-props
:open:
We verify some of the fundamental SVD properties using the built-in `svd` function.

```{index} ! MATLAB; svd
```

```{code-cell}
A = vander(1:5);
A = A(:, 1:4)
```

```{code-cell}
[U, S, V] = svd(A);
disp(sprintf("U is %d by %d. S is %d by %d. V is %d by %d.\n", size(U), size(S), size(V)))
```

We verify the orthogonality of the singular vectors as follows:

```{code-cell}
norm(U' * U - eye(5))
norm(V' * V - eye(4))
```

Here is verification of the connections between the singular values, norm, and condition number.

```{code-cell}
s = diag(S);
norm_A = norm(A)
sigma_max = s(1)
```

```{code-cell}
cond_A = cond(A)
sigma_ratio = s(1) / s(end)
```
``````

### 7.4 @section-matrixanaly-symm-eig

(demo-symm-eig-rayleigh-matlab)=
``````{dropdown} @demo-symm-eig-rayleigh
:open:
We will use a symmetric matrix with a known EVD and eigenvalues equal to the integers from 1 to 20.

```{code-cell}
n = 20;
lambda = 1:n;
D = diag(lambda);
[V, ~] = qr(randn(n, n));    % get a random orthogonal V
A = V * D * V';
```

The Rayleigh quotient is a scalar-valued function of a vector.

```{code-cell}
R = @(x) (x' * A * x) / (x' * x);
```

The Rayleigh quotient evaluated at an eigenvector gives the corresponding eigenvalue.

```{code-cell}
format long
R(V(:, 7))
```

If the input to he Rayleigh quotient is within a small $\delta$ of an eigenvector, its output is within $O(\delta^2)$ of the corresponding eigenvalue. In this experiment, we observe that each additional digit of accuracy in an approximate eigenvector gives two more digits to the eigenvalue estimate coming from the Rayleigh quotient.

```{code-cell}
delta = 1 ./ 10 .^ (1:5)';
dif = zeros(size(delta));
for k = 1:length(delta)
    e = randn(n, 1);
    e = delta(k) * e / norm(e);
    x = V(:, 6) + e;
    dif(k) = R(x) - lambda(6);
end
table(delta, dif, variablenames=["perturbation size", "R(x) - lambda"])
```
``````

### 7.5 @section-matrixanaly-dimreduce

(demo-dimreduce-hello-matlab)=
``````{dropdown} @demo-dimreduce-hello
:open:
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
``````

(demo-dimreduce-voting-matlab)=
``````{dropdown} @demo-dimreduce-voting
:open:
This matrix describes the votes on bills in the 111th session of the United States Senate. (The data set was obtained from [https://voteview.com].) Each row is one senator, and each column is a vote item.

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
``````