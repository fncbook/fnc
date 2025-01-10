---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---
# Chapter 8

## Functions

(function-poweriter-matlab)=
``````{dropdown} Power iteration
:open:
```{literalinclude} FNC-matlab/poweriter.m
:language: matlab
:linenos: true
```
``````

(function-inviter-matlab)=
``````{dropdown} Inverse iteration
:open:
```{literalinclude} FNC-matlab/inviter.m
:language: matlab
:linenos: true
```
``````

(function-arnoldi-matlab)=
``````{dropdown} Arnoldi iteration
:open:
```{literalinclude} FNC-matlab/arnoldi.m
:language: matlab
:linenos: true
```
``````

(function-gmres-matlab)=
``````{dropdown} GMRES
:open:
```{literalinclude} FNC-matlab/arngmres.m
:language: matlab
:linenos: true
```
``````

## Examples

```{code-cell}
:tags: [remove-cell]
cd  /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init
pwd
```

### 8.1 @section-krylov-structure

(demo-structure-sparse-matlab)=
``````{dropdown} @demo-structure-sparse
Here we load the adjacency matrix of a graph with 2790 nodes. Each node is a web page referring to Roswell, NM, and the edges represent links between web pages. (Credit goes to Panayiotis Tsaparas and the University of Toronto for making this data public.)

```{code-cell}
load roswelladj
a = whos('A')
```

```{index} ! MATLAB; nnz
```
We may define the density of $\mathbf{A}$ as the number of nonzeros divided by the total number of entries.
```{tip}
:class: dropdown
Use `nnz` to count the number of nonzeros in a sparse matrix.
```

```{code-cell}
sz = size(A);  n = sz(1);
density = nnz(A) / prod(sz)
```

The computer memory consumed by any variable can be discovered using `whos`. We can use it to compare the space needed for the sparse representation to its dense counterpart, that is, the space needed to store all the elements, whether zero or not.

```{code-cell}
F = full(A);
f = whos('F');
storage_ratio = f.bytes / a.bytes
```

Matrix-vector products are also much faster using the sparse form because operations with structural zeros are skipped.

```{code-cell}
x = randn(n,1);
tic, for i = 1:200, A*x; end
sparse_time = toc
```

```{code-cell}
tic, for i = 1:200, F*x; end
dense_time = toc
```

However, the sparse storage format in MATLAB is column-oriented.  Operations on rows may take a lot longer than similar ones on columns.

```{code-cell}
v = A(:, 1000);
tic, for i = 1:n, A(:, i) = v; end
column_time = toc
r = v';
tic, for i = 1:n, A(i, :) = r; end
row_time = toc
``````

(demo-structure-fill-matlab)=
``````{dropdown} @demo-structure-fill

Here is the adjacency matrix of a graph representing a small-world network, featuring connections to neighbors and a small number of distant contacts.

```{code-cell}
load smallworld.mat
G = graph(A);
plot(G, nodecolor='r')
```

Because each node connects to relatively few others, the adjacency matrix is quite sparse.

```{code-cell}
spy(A)
```

By {numref}`Theorem {number} <theorem-insight-adjmat>`, the entries of $\mathbf{A}^k$ give the number of walks of length $k$ between pairs of nodes, as with "*k* degrees of separation" within a social network. As $k$ grows, the density of $\mathbf{A}^k$ also grows.

```{code-cell}
clf
tiledlayout(2, 2)
for k = [2, 3, 4, 6]
    nexttile
    spy(A^k)
    title(sprintf("A^{%d}", k))
end
```
``````

(demo-structure-sparseband-matlab)=
``````{dropdown} @demo-structure-sparseband 

```{index} ! MATLAB; spdiags
```

The `spdiags` function creates a sparse matrix given its diagonal elements. The main or central diagonal is numbered zero, above and to the right of that is positive, and below and to the left is negative.

```{code-cell}
n = 50;
n = 50;
% Put constant values on 3 diagonals
A = spdiags([n, 1, 0.1], [-3, 0, 5], n, n);
% Put other values on 1st superdiagonal
A = spdiags(-(0:n-1)', 1, A);
full(A(1:7, 1:7))
```

```{index} ! MATLAB; sparse
```

Without pivoting, the LU factors have the same lower and upper bandwidth as the original matrix.
```{tip}
:class: dropdown
The `sparse` function converts any matrix to sparse form. But it's usually better to construct a sparse matrix directly, as the standard form might not fit in memory.
```



```{code-cell}
[L, U] = lufact(A);
clf
subplot(1, 2, 1), spy(L), title('L')
subplot(1, 2, 2), spy(U), title(('U'));
```

However, if we introduce row pivoting, bandedness may be expanded or destroyed.

```{code-cell}
[L, U, p] = plufact(A);
subplot(1, 2, 1), spy(L(p, :)), title('L')
subplot(1, 2, 2), spy(U), title(('U'));
```
``````

(demo-structure-linalg-matlab)=

``````{dropdown} @demo-structure-linalg

The following generates a random sparse matrix with prescribed eigenvalues.

```{code-cell}
n = 4000;
density = 4e-4;
lambda = 1 ./ (1:n);
A = sprandsym(n, density, lambda);
clf,  spy(A)
title('Sparse symmetric matrix')  
```

```{index} ! MATLAB; eigs
```

The `eigs` function finds a small number eigenvalues meeting some criterion. First, we ask for the 5 of largest (complex) magnitude.

```{code-cell}
[V, D] = eigs(A, 5);    % largest magnitude
1 ./ diag(D)            % should be 1, 2, 3, 4, 5
```

Now we find the 4 closest to the value 0.03 in the complex plane.

```{code-cell}
[V, D] = eigs(A, 4, 0.03);    % closest to 0.03
diag(D)
```

```{index} MATLAB; \\
```

The time needed to solve a sparse linear system is not easy to predict unless you have some more information about the matrix. But it will typically be orders of magnitude faster than the dense version of the same problem.

```{code-cell}
x = 1 ./ (1:n)';  
b = A * x;
tic, sparse_err = norm(x - A\b), sparse_time = toc
```

```{code-cell}
F = full(A);
tic, dense_err = norm(x - F\b), dense_time = toc
```
``````

### 8.2 @section-krylov-power

(demo-power-one-matlab)=
``````{dropdown} @demo-power-one
Here we choose a magic 5×5 matrix and a random 5-vector.

```{code-cell}
A = magic(5) / 65;
x = randn(5, 1);
```

Applying matrix-vector multiplication once doesn't do anything recognizable.

```{code-cell}
y = A * x
```

Repeating the multiplication still doesn't do anything obvious.

```{code-cell}
z = A * y
```

But if we keep repeating the matrix-vector multiplication, something remarkable happens: $\mathbf{A} \mathbf{x} \approx \mathbf{x}$.

```{code-cell}
for j = 1:8
    x = A * x;
end
[x, A * x]
```

This phenomenon seems to occur regardless of the starting vector.

```{code-cell}
x = randn(5, 1);
for j = 1:8
    x = A * x;
end
[x, A * x]
```
``````

(demo-power-iter-matlab)=
``````{dropdown} @demo-power-iter
We will experiment with the power iteration on a 5×5 matrix with prescribed eigenvalues and dominant eigenvalue at 1.

```{code-cell}
ev = [1, -0.75, 0.6, -0.4, 0];
% Make a triangular matrix with eigenvalues on the diagonal.
A = triu(ones(5, 5), 1) + diag(ev);
```

We run the power iteration 60 times. The best estimate of the dominant eigenvalue is the last entry of the first output.

```{code-cell}
[beta, x] = poweriter(A, 60);
format long
beta(1:12)
```

We check for linear convergence using a log-linear plot of the error.

```{code-cell}
err = 1 - beta;
clf,  semilogy(abs(err), '.-')
title('Convergence of power iteration')
xlabel('k'),  ylabel(('|\lambda_1 - \beta_k|'));
```

The asymptotic trend seems to be a straight line, consistent with linear convergence. To estimate the convergence rate, we look at the ratio of two consecutive errors in the linear part of the convergence curve. The ratio of the first two eigenvalues should match the observed rate.

```{code-cell}
theory = ev(2) / ev(1)
observed = err(40) / err(39)
```

Note that the error is supposed to change sign on each iteration. The effect of these alternating signs is that estimates oscillate around the exact value.

```{code-cell}
beta(26:29)
```

In practical situations, we don't know the exact eigenvalue that the algorithm is supposed to find. In that case we would base errors on the final $\beta$ that was found, as in the following plot.

```{code-cell}
err = beta(end) - beta(1:end-1);
semilogy(abs(err), '.-')
title('Convergence of power iteration')
xlabel('k'),  ylabel(('|\beta_{60} - \beta_k|'));
```

The results are very similar until the last few iterations, when the limited accuracy of the reference value begins to show. That is, while it is a good estimate of $\lambda_1$, it is less good as an estimate of the error in nearby estimates.
``````


### 8.3 @section-krylov-inviter

(demo-inviter-conv-matlab)=
``````{dropdown} @demo-inviter-conv
We set up a $5\times 5$ triangular matrix with prescribed eigenvalues on its diagonal.

```{code-cell}
ev = [1, -0.75, 0.6, -0.4, 0];
A = triu(ones(5, 5), 1) + diag(ev);
```

We run inverse iteration with the shift $s=0.7$. The result should converge to the eigenvalue closest to 0.7, which we know to be 0.6 here.

```{code-cell}
s = 0.7;
[beta, x] = inviter(A, s, 30);
format short
beta(1:10)
```

The convergence is again linear.

```{code-cell}
err = abs(0.6 - beta);
semilogy(abs(err),'.-')
title('Convergence of inverse iteration')
xlabel('k'), ylabel(('|\lambda_j - \beta_k|'));
```

Let's reorder the eigenvalues to enforce {eq}`shiftorder`.
```{tip}
:class: dropdown
The second output of `sort` returns the index permutation needed to sort the given vector.
```

```{code-cell}
[~, idx] = sort(abs(ev - s));
ev = ev(idx)
```

Now it is easy to compare the theoretical and observed linear convergence rates.

```{code-cell}
theoretical_rate = (ev(1) - s) / (ev(2) - s)
observed_rate = err(26) / err(25)
```
``````

(demo-inviter-accel-matlab)=
``````{dropdown} @demo-inviter-accel
```{code-cell}
ev = [1, -0.75, 0.6, -0.4, 0];
A = triu(ones(5, 5), 1) + diag(ev);
```

We begin with a shift $s=0.7$, which is closest to the eigenvalue 0.6.

```{code-cell}
s = 0.7;
x = ones(5, 1);
y = (A - s * eye(5)) \ x;
beta = x(1) / y(1) + s
```

Note that the result is not yet any closer to the targeted 0.6. But we proceed (without being too picky about normalization here).

```{code-cell}
s = beta;
x = y / y(1);
y = (A - s * eye(5)) \ x;
beta = x(1) / y(1) + s
```

Still not much apparent progress. However, in just a few more iterations the results are dramatically better.

```{code-cell}
format long
for k = 1:4
    s = beta;
    x = y / y(1);
    y = (A - s * eye(5)) \ x;
    beta = x(1) / y(1) + s
end
```
``````

### 8.4 @section-krylov-subspace

(demo-subspace-unstable-matlab)=
``````{dropdown} @demo-subspace-unstable
First we define a triangular matrix with known eigenvalues, and a random vector $b$.

```{code-cell}
lambda = 10 + (1:100);
A = diag(lambda) + triu(rand(100), 1); 
b = rand(100, 1);
```

Next we build up the first ten Krylov matrices iteratively, using renormalization after each matrix-vector multiplication.

```{code-cell}
Km = b;
for m = 1:29      
    v = A * Km(:, m);
    Km(:, m+1) = v / norm(v);
end
```

Now we solve least-squares problems for Krylov matrices of increasing dimension, recording the residual in each case.

```{code-cell}
warning off  
resid = zeros(30, 1);
for m = 1:30  
    z = (A * Km(:, 1:m)) \ b;
    x = Km(:, 1:m) * z;
    resid(m) = norm(b - A * x);
end
```

The linear system approximations show smooth linear convergence at first, but the convergence stagnates after only a few digits have been found.

```{code-cell}
clf
semilogy(resid, '.-')
xlabel('m'),  ylabel('|| b-Ax_m ||')
set(gca,'ytick',10.^(-6:2:0))
axis tight, title('Residual for linear systems') 
```
``````

(demo-subspace-arnoldi-matlab)=
``````{dropdown} @demo-subspace-arnoldi
We illustrate a few steps of the Arnoldi iteration for a small matrix.

```{code-cell}
A = randi(9, [6, 6])
```

The seed vector we choose here determines the first member of the orthonormal basis.

```{code-cell}
u = randn(6, 1);
Q = u / norm(u);
```

Multiplication by $\mathbf{A}$ gives us a new vector in $\mathcal{K}_2$.

```{code-cell}
Aq = A * Q(:, 1);
```

We subtract off its projection in the previous direction. The remainder is rescaled to give us the next orthonormal column.

```{code-cell}
v = Aq - dot(Q(:, 1), Aq) * Q(:, 1);
Q = [Q, v / norm(v)];
```

On the next pass, we have to subtract off the projections in two previous directions.

```{code-cell}
Aq = A * Q(:, 2);
v = Aq - dot(Q(:, 1), Aq) * Q(:, 1) - dot(Q(:, 2), Aq) * Q(:, 2);
Q = [Q, v / norm(v)];
```

At every step, $\mathbf{Q}_m$ is an ONC matrix.

```{code-cell}
format
norm(Q' * Q - eye(3))
```

And $\mathbf{Q}_m$ spans the same space as the three-dimensional Krylov matrix.

```{code-cell}
K = [u, A * u, A * A * u];
rank([Q, K])
```
``````

### 8.5 @section-krylov-gmres

(demo-gmres-intro-matlab)=
``````{dropdown} @demo-gmres-intro
We define a triangular matrix with known eigenvalues and a random vector $\mathbf{b}$.

```{code-cell}
lambda = 10 + (1:100);
A = diag(lambda) + triu(rand(100), 1); 
b = rand(100, 1);
```

Instead of building the Krylov matrices, we use the Arnoldi iteration to generate equivalent orthonormal vectors.

```{code-cell}
[Q, H] = arnoldi(A, b, 60);
```

The Arnoldi bases are used to solve the least-squares problems defining the GMRES iterates.

```{code-cell}
resid = norm(b);
for m = 1:60
    s = [norm(b); zeros(m, 1)];
    z = H(1:m+1, 1:m) \ s;
    x = Q(:, 1:m) * z;
    resid = [resid, norm(b - A * x)];
end
```

The approximations converge smoothly, practically all the way to machine epsilon.

```{code-cell}
clf
semilogy(resid,'.-')
xlabel('m'),  ylabel('|| b - Ax_m ||')
axis tight, title(('Residual for GMRES'));
```
``````

(demo-gmres-restart-matlab)=
``````{dropdown} @demo-gmres-restart
The following experiments are based on a matrix resulting from discretization of a partial differential equation.

```{code-cell}
d = 50;
A = d^2 * gallery('poisson', d);
n = size(A, 1)
b = ones(n, 1);
clf,  spy(A)
```

```{index} ! MATLAB; gmres
```

We compare unrestarted GMRES with three different thresholds for restarting. Here we are using the built-in `gmres`, since our simple implementation does not offer restarting.

```{code-cell}
clf
restart = [120, 20, 40, 60];
for j = 1:4
    [~,~,~,~,rv] = gmres(A, b, restart(j), 1e-9,120 / restart(j));
    semilogy(0:length(rv) - 1, rv),  hold on
end
title('Convergence of restarted GMRES')
xlabel('m'),  ylabel('residual norm')
legend('no restart','every 20','every 40','every 60','location','southwest');
```

The "pure" GMRES curve is the lowest one. All of the other curves agree with it until the first restart. Decreasing the restart value makes the convergence per iteration generally worse, but the time required per iteration smaller as well.

``````
### 8.6 @section-krylov-minrescg

(demo-minrescg-indefinite-matlab)=
``````{dropdown} @demo-minrescg-indefinite
The following matrix is indefinite.

```{code-cell}
A = (11 / pi)^2 * gallery('poisson', 10);
A = A - 20 * eye(100);
lambda = eig(full(A));
isneg = lambda < 0;
disp(sprintf("%d negative and %d positive eigenvalues", sum(isneg), sum(~isneg)))
```

We can compute the relevant quantities from {numref}`Theorem {number} <theorem-minrescg-indefinite>`.

```{code-cell}
m = min(-lambda(isneg));
M = max(-lambda(isneg));
kappa_minus = M / m;
m = min(lambda(~isneg));
M = max(lambda(~isneg));
kappa_plus = M / m;
S = sqrt(kappa_minus * kappa_plus);
rho = sqrt((S - 1) / (S + 1));
fprintf("convergence rate: %.3f", rho)
```

Because the iteration number $m$ is halved in {eq}`minres-conv`, the rate constant of linear convergence is the square root of this number, which makes it even closer to 1.

Now we apply MINRES to a linear system with this matrix, and compare the observed convergence to the upper bound from the theorem.

```{index} ! MATLAB; minres
```

```{code-cell}
b = rand(100, 1);
[xMR, ~,~ , ~, residMR] = minres(A, b, 1e-10, 100);
relres = residMR / norm(b);
m = 0:length(relres) - 1;
clf,  semilogy(m, relres, '.-')
hold on
semilogy(m, rho .^ m, 'k--')
xlabel('m'),  ylabel('relative residual') 
title('Convergence of MINRES') 
legend('MINRES', 'upper bound', 'location', 'southwest');
```

The upper bound turns out to be pessimistic here, especially in the later iterations. However, you can see a slow linear phase in the convergence that tracks the bound closely.
``````

(demo-minrescg-converge-matlab)=
``````{dropdown} @demo-minrescg-converge
We will compare MINRES and CG on some quasi-random SPD problems.  The first matrix has a condition number of 100.

```{code-cell}
n = 5000;
density = 0.001;
A = sprandsym(n, density, 1e-2, 2);
```

We generate a system with a known solution.

```{code-cell}
x = (1:n)' / n;
b = A * x;
```

```{index} ! MATLAB; pcg
```

Now we apply both methods and compare the convergence of the system residuals.

```{code-cell}
[xMR, ~, ~, ~, residMR] = minres(A, b, 1e-7, 100);
[xCG, ~, ~, ~, residCG] = pcg(A, b, 1e-7, 100);
M = length(residMR) - 1;
clf,  semilogy(0:M, residMR / norm(b), '.-')
M = length(residCG) - 1;
hold on,  semilogy(0:M, residCG / norm(b), '.-')
title('Convergence of MINRES and CG')
xlabel('Krylov dimension m')
ylabel('||r_m|| / ||b||')
legend('MINRES', 'CG');
```

There is little difference between the two methods here. Next, we increase the condition number of the matrix by a factor of 25. The rule of thumb predicts that the number of iterations required should increase by a factor of about 5.

```{code-cell}
A = sprandsym(n, density, 1e-2 / 25, 2);
b = A * x;
```

```{code-cell}
:tags: [hide-input]

[xMR, ~, ~, ~, residMR] = minres(A, b, 1e-7, 400);
[xCG, ~, ~, ~, residCG] = pcg(A, b, 1e-7, 400);
M = length(residMR) - 1;
clf,  semilogy(0:M, residMR / norm(b), '.-')
M = length(residCG) - 1;
hold on,  semilogy(0:M, residCG / norm(b), '.-')
title('Convergence of MINRES and CG')
xlabel('Krylov dimension m')
ylabel('||r_m|| / ||b||')
legend('MINRES', 'CG');
```

Both methods have an early superlinear phase that allow them to finish slightly sooner than the factor of 5 predicted: {numref}`Theorem {number} <theorem-minrescg-converge>` is an upper bound, not necessarily an approximation. Both methods ultimately achieve the same reduction in the residual; MINRES stops earlier, but with a slightly larger error.

``````

### 8.7 @section-krylov-matrixfree

(demo-matrixfree-blur-matlab)=
``````{dropdown} @demo-matrixfree-blur
We use a readily available test image.

```{code-cell}
load mandrill
[m, n] = size(X);
clf
imshow(X, [0, 255])
title('Original image')    % ignore this 
```

We define the one-dimensional tridiagonal blurring matrices.

```{code-cell}
v = [1/4, 1/2, 1/4];
B = spdiags(v, -1:1, m, m);
C = spdiags(v, -1:1, n, n);
```

Finally, we show the results of using $k=12$ repetitions of the blur in each direction.

```{code-cell}
blur = @(X) B^12 * X * C^12;
imshow(blur(X), [0, 255])
title(('Blurred image'));
```
``````

(demo-matrixfree-deblur-matlab)=
``````{dropdown} @demo-matrixfree-deblur
We repeat the earlier process to blur an original image $\mathbf{X}$ to get $\mathbf{Z}$.

```{code-cell}
:tags: hide-input
load mandrill
[m, n] = size(X);
v = [1/4, 1/2, 1/4];
B = spdiags(v, -1:1, m, m);
C = spdiags(v, -1:1, n, n);
blur = @(X) B^12 * X * C^12;
```

```{code-cell}
Z = blur(X);
clf,  imshow(Z, [0, 255])
title(("Blurred image"));
```

Now we imagine that $\mathbf{X}$ is unknown and that we want to recover it from $\mathbf{Z}$. We first need functions that translate between vector and matrix representations.

```{code-cell}
vec = @(X) reshape(X,m*n,1);
unvec = @(x) reshape(x,m,n);
T = @(x) vec( blur(unvec(x)) );
```
The blurring operators are symmetric, so we apply `minres` to the composite blurring transformation `T`.

```{code-cell}
y = gmres(T, vec(Z), 50, 1e-5);
Y = unvec(y);

subplot(121)
imshow(X, [0, 255])
title("Original")
subplot(122)
imshow(Y, [0, 255])
title(("Deblurred"));
```
``````

### 8.8 @section-krylov-precond

(demo-precond-diagonal-matlab)=
``````{dropdown} @demo-precond-diagonal
Here is an SPD matrix that arises from solving partial differential equations.

```{code-cell}
A = gallery("wathen", 60, 60);
n = size(A, 1);
clf,  spy(A)
```

There is an easy way to use the diagonal elements of $\mathbf{A}$, or any other vector, as a diagonal preconditioner.

```{code-cell}
M = spdiags(diag(A), 0, n, n);
```

We now compare MINRES with and without the preconditioner.

```{code-cell}
:tags: hide-input

b = ones(n, 1);
[x, ~, ~, ~, resid_plain] = minres(A, b, 1e-10, 400);
clf,  semilogy(resid_plain)
xlabel('iteration number'), ylabel('residual norm')
title('Unpreconditioned MINRES')

[x, ~, ~, ~, resid_prec] = minres(A, b, 1e-10, 400, M);
hold on,  semilogy(resid_prec)
title('Precondtioned MINRES')
legend('no prec.', 'with prec.');
```

The diagonal preconditioner cut down substantially on the number of iterations. The effect on the total time is less dramatic, but this is not a large version of the problem.

``````

(demo-precond-gmres-matlab)=
``````{dropdown} @demo-precond-gmres
Here is a random nonsymmetric matrix.

```{code-cell}
n = 8000;
A = speye(n) + sprand(n, n, 0.00035);
```

Without a preconditioner, restarted GMRES makes slow progress.

```{code-cell}
b = rand(n, 1);
[x, ~, ~, ~, resid_plain] = gmres(A, b, 50, 1e-10, 3);  % restart at 50
format short e
resid_plain(1:30:end)
```

```{index} ! MATLAB; ilu
```

This version of incomplete LU factorization simply prohibits fill-in for the factors, freezing the sparsity pattern of the approximate factors to match the original matrix.

```{code-cell}
[L, U] = ilu(A);
clf
subplot(121), spy(L)
title('L')
subplot(122), spy(U)
title('U')
disp(sprintf("There are %d nonzeros in A", nnz(A)))
```

It does _not_ produce a true factorization of $\mathbf{A}$.

```{code-cell}
norm( full(A - L * U) )  
```

The actual preconditioning matrix is $\mathbf{M}=\mathbf{L}\mathbf{U}$. However, the `gmres` function allows setting the preconditioner by giving the factors independently.

```{code-cell}
[x, ~, ~, ~, resid_prec] = gmres(A, b, [], 1e-10, 300, L, U);
```

The preconditioner makes a significant difference in the number of iterations needed.

```{code-cell}
clf, semilogy(resid_plain)
hold on, semilogy(resid_prec)
xlabel('iteration number'), ylabel('residual norm')
title('Precondtioned GMRES ')
legend('no preconditioner', 'with preconditioner');
```
```