---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
# Chapter 8

## Functions

(function-poweriter-python)=
``````{dropdown} Power iteration
:open:
```{literalinclude} fncbook/fncbook/chapter08.py
:filename: poweriter.py
:start-at: def poweriter
:end-at: return gamma, x
:language: python
:linenos: true
```
``````

(function-inviter-python)=
``````{dropdown} Inverse iteration
:open:
```{literalinclude} fncbook/fncbook/chapter08.py
:filename: inviter.py
:start-at: def inviter
:end-at: return gamma, x
:language: python
:linenos: true
```
``````

(function-arnoldi-python)=
``````{dropdown} Arnoldi iteration
:open:
```{literalinclude} fncbook/fncbook/chapter08.py
:filename: arnoldi.py
:start-at: def arnoldi
:end-at: return Q, H
:language: python
:linenos: true
```

```{admonition} About the code
The loop starting at line 17 does not exactly implement {eq}`arnoldiip` and {eq}`arnoldigs`. The reason is numerical stability. Though the described and implemented versions are mathematically equivalent in exact arithmetic (see [Exercise 6](#problem-subspace-modifiedgs)), the approach in {numref}`Function {number} <function-arnoldi>` is more stable.
```
``````

(function-gmres-python)=
``````{dropdown} GMRES
:open:
```{literalinclude} fncbook/fncbook/chapter08.py
:filename: gmres.py
:start-at: def gmres
:end-at: return x, residual
:language: python
:linenos: true
```
``````

## Examples

```{code-cell} ipython3
:tags: remove-cell
exec(open("FNC_init.py").read())
```

### 8.1 @section-krylov-structure

(demo-structure-sparse-python)=
``````{dropdown} @demo-structure-sparse
:open:
```{tip}
Functions to work with sparse matrices are found in the `scipy.sparse` module.
```

Here we load the adjacency matrix of a graph with 2790 nodes. Each node is a web page referring to Roswell, NM, and the edges represent links between web pages. (Credit goes to Panayiotis Tsaparas and the University of Toronto for making this data public.)

```{code-cell}
import scipy.sparse as sp
from scipy.io import loadmat

vars = loadmat("roswelladj.mat")    # get from the book's website
A = sp.csr_matrix(vars["A"])
```

We may define the density of $\mathbf{A}$ as the number of nonzeros divided by the total number of entries.

```{index} ! Python; nnz
```

```{code-cell}
m, n = A.shape
print(f"density is {A.nnz / (m * n):.3%}")
```

We can compare the storage space needed for the sparse $\mathbf{A}$ with the space needed for its dense / full counterpart.


```{code-cell}
F = A.todense()
print(f"{A.data.nbytes/1e6:.3f} MB for sparse form, {F.nbytes/1e6:.3f} MB for dense form")
```

Matrix-vector products are also much faster using the sparse form because operations with structural zeros are skipped.

```{code-cell}
from timeit import default_timer as timer
x = random.randn(n)
start = timer()
for i in range(1000):
    A @ x
print(f"sparse time: {timer() - start:.4g} sec")
```

```{code-cell}
start = timer()
for i in range(1000):
    F @ x
print(f"dense time: {timer() - start:.4g} sec")
```
``````

(demo-structure-fill-python)=
``````{dropdown} @demo-structure-fill
:open:

Here is the adjacency matrix of a graph representing a small-world network, featuring connections to neighbors and a small number of distant contacts.

```{code-cell}
import networkx as nx
wsg = nx.watts_strogatz_graph(200, 4, 0.02)
```

Because each node connects to relatively few others, the adjacency matrix is quite sparse.

```{code-cell}
A = nx.adjacency_matrix(wsg)
spy(A)
title("Adjacency matrix $A$");
```

By {numref}`Theorem {number} <theorem-insight-adjmat>`, the entries of $\mathbf{A}^k$ give the number of walks of length $k$ between pairs of nodes, as with "*k* degrees of separation" within a social network. As $k$ grows, the density of $\mathbf{A}^k$ also grows.
```{tip}
:class: dropdown
While `A**6` is valid syntax here, it means elementwise power, not matrix power. 
```

```{index} ! Python; matrix_power
```

```{code-cell}
from scipy.sparse.linalg import matrix_power
spy(matrix_power(A, 6))
title(("$A^6$"));
```
``````

(demo-structure-sparseband-python)=
``````{dropdown} @demo-structure-sparseband
:open:
```{index} ! Julia; spdiagm
```

The `scipi.sparse.diags` function creates a sparse matrix given its diagonal elements and the diagonal indexes to put them on. The main or central diagonal is numbered zero, above and to the right of that is positive, and below and to the left is negative.

```{code-cell}
n = 50
data = [n * ones(n-3), ones(n), linspace(-1, 1-n, n-1)]
offsets = [-3, 0, 1]    # 3rd below, main, 1st above
A = sp.diags(data, offsets, format="lil")
print(A[:7, :7].todense())
```

Without pivoting, the LU factors have the same lower and upper bandwidth as the original matrix.

```{code-cell}
L, U = FNC.lufact(A.todense())
subplot(1, 2, 1), spy(L)
subplot(1, 2, 2), spy(U);
```

However, if we introduce row pivoting, bandedness may be expanded or destroyed.

```{code-cell}
L, U, p = FNC.plufact(A.todense())
subplot(1, 2, 1), spy(L[p, :])
subplot(1, 2, 2), spy(U)
```
``````

(demo-structure-linalg-python)=

``````{dropdown} @demo-structure-linalg
:open:

The following generates a random sparse matrix with prescribed eigenvalues.

```{code-cell}
n = 4000
density = 4e-4
ev = 1 / arange(1, n + 1)
A = FNC.sprandsym(n, density, eigvals=ev)
print(f"density is {A.nnz / prod(A.shape):.3%}")
```

```{index} ! Python; eigs
```

The `eigs` function finds a small number eigenvalues meeting some criterion. First, we ask for the 5 of largest (complex) magnitude using `which="LM"`.

```{code-cell}
from scipy.sparse.linalg import eigs
ev, V = eigs(A, k=5, which="LM")    # largest magnitude
print(1 / ev)
```

Now we find the 4 closest to the value 1 in the complex plane, via `sigma=1`.

```{code-cell}
from scipy.sparse.linalg import eigs
ev, V = eigs(A, k=4, sigma=0.03)    # closest to sigma
print(ev)
```

The time needed to solve a sparse linear system is not easy to predict unless you have some more information about the matrix. But it will typically be orders of magnitude faster than the dense version of the same problem.

```{code-cell}
from scipy.sparse.linalg import spsolve
x = 1 / arange(1, n + 1)
b = A @ x
start = timer()
xx = spsolve(A, b)
print(f"sparse time: {timer() - start:.3g} sec")
print(f"residual: {norm(b - A @ xx, 2):.1e}")
```

```{code-cell}
from numpy.linalg import solve
F = A.todense()
start = timer()
xx = solve(F, b)
print(f"dense time: {timer() - start:.3g} sec")
print(f"residual: {norm(b - A @ xx, 2):.1e}")
```
``````

### 8.2 @section-krylov-power

(demo-power-one-python)=
``````{dropdown} @demo-power-one
:open:
Here we choose a random 5×5 matrix and a random 5-vector.

```{code-cell}
A = random.choice(range(10), (5, 5))
A = A / sum(A, 0)
x = random.randn(5)
print(x)
```

Applying matrix-vector multiplication once doesn't do anything recognizable.

```{code-cell}
y = A @ x
print(y)
```

Repeating the multiplication still doesn't do anything obvious.

```{code-cell}
z = A @ y
print(z)
```

But if we keep repeating the matrix-vector multiplication, something remarkable happens: $\mathbf{A} \mathbf{x} \approx \mathbf{x}$.

```{code-cell}
x = random.randn(5)
for j in range(6):
    x = A @ x
print(x)
print(A @ x)
```

This phenomenon is unlikely to be a coincidence!

``````

(demo-power-iter-python)=
``````{dropdown} @demo-power-iter
:open:
We will experiment with the power iteration on a 5×5 matrix with prescribed eigenvalues and dominant eigenvalue at 1.

```{code-cell}
ev = [1, -0.75, 0.6, -0.4, 0]
A = triu(ones([5, 5]), 1) + diag(ev)    # triangular matrix, eigs on diagonal
```

We run the power iteration 60 times. The first output should be a sequence of estimates converging to the dominant eigenvalue—which, in this case, we set up to be 1.

```{code-cell}
beta, x = FNC.poweriter(A, 60)
print(beta)
```

We check for linear convergence using a log-linear plot of the error.

```{code-cell}
err = 1 - beta
semilogy(arange(60), abs(err), "-o")
ylim(1e-10, 1)
xlabel("$k$")
ylabel("$|\\lambda_1 - \\beta_k|$")
title("Convergence of power iteration");
```

The asymptotic trend seems to be a straight line, consistent with linear convergence. To estimate the convergence rate, we look at the ratio of two consecutive errors in the linear part of the convergence curve. The ratio of the first two eigenvalues should match the observed rate.

```{code-cell}
print(f"theory: {ev[1] / ev[0]:.5f}")
print(f"observed: {err[40] / err[39]:.5f}")
```

Note that the error is supposed to change sign on each iteration. The effect of these alternating signs is that estimates oscillate around the exact value.

```{code-cell}
print(beta[26:30])
```

In practical situations, we don't know the exact eigenvalue that the algorithm is supposed to find. In that case we would base errors on the final $\beta$ that was found, as in the following plot.

```{code-cell}
err = beta[-1] - beta
semilogy(arange(60), abs(err), "-o")
ylim(1e-10, 1), xlabel("$k$")
ylabel("$|\\lambda_1 - \\beta_k|$")
title("Convergence of power iteration");
```

The results are very similar until the last few iterations, when the limited accuracy of the reference value begins to show. That is, while it is a good estimate of $\lambda_1$, it is less good as an estimate of the error in nearby estimates.

``````

### 8.3 @section-krylov-inviter

(demo-inviter-conv-python)=
``````{dropdown} @demo-inviter-conv
:open:
We set up a $5\times 5$ triangular matrix with prescribed eigenvalues on its diagonal.

```{code-cell}
ev = array([1, -0.75, 0.6, -0.4, 0])
A = triu(ones([5, 5]), 1) + diag(ev)    # triangular matrix, eigs on diagonal
```

We run inverse iteration with the shift $s=0.7$. Convergence should be to the eigenvalue closest to the shift, which we know to be $0.6$ here.

```{code-cell}
beta, x = FNC.inviter(A, 0.7, 30)
print(beta)
```

As expected, the eigenvalue that was found is the one closest to 0.7. The convergence is again linear.

```{code-cell}
err = beta[-1] - beta    # last estimate is our best
semilogy(arange(30), abs(err), "-o")
ylim(1e-16, 1)
xlabel("$k$"),  ylabel("$|\\lambda_3 - \\beta_k|$")
title(("Convergence of inverse iteration"));
```

```{index} ! Python; argsort
```

Let's reorder the eigenvalues to enforce {eq}`shiftorder`.
```{tip}
:class: dropdown
The `argsort` function returns the index permutation needed to sort the given vector, rather than the sorted vector itself.
```

```{code-cell}
ev = ev[argsort(abs(ev - 0.7))]
print(ev)
```

Now it is easy to compare the theoretical and observed linear convergence rates.

```{code-cell}
print(f"theory: {(ev[0] - 0.7) / (ev[1] - 0.7):.5f}")
print(f"observed: {err[21] / err[20]:.5f}")
```
``````

(demo-inviter-accel-python)=
``````{dropdown} @demo-inviter-accel
:open:
```{code-cell}
ev = array([1, -0.75, 0.6, -0.4, 0])
A = triu(ones([5, 5]), 1) + diag(ev)    # triangular matrix, eigs on diagonal
```

We begin with a shift $s=0.7$, which is closest to the eigenvalue 0.6.

```{code-cell}
from numpy.linalg import solve
s = 0.7
x = ones(5)
y = solve(A - s * eye(5), x)
beta = x[0] / y[0] + s
print(f"latest estimate: {beta:.8f}")
```

Note that the result is not yet any closer to the targeted 0.6. But we proceed (without being too picky about normalization here).

```{code-cell}
s = beta
x = y / y[0]
y = solve(A - s * eye(5), x)
beta = x[0] / y[0] + s
print(f"latest estimate: {beta:.8f}")
```

Still not much apparent progress. However, in just a few more iterations the results are dramatically better.

```{code-cell}
for k in range(4):
    s = beta
    x = y / y[0]
    y = solve(A - s * eye(5), x)
    beta = x[0] / y[0] + s
    print(f"latest estimate: {beta:.12f}")
```
``````

### 8.4 @section-krylov-subspace

(demo-subspace-unstable-python)=
``````{dropdown} @demo-subspace-unstable
:open:
First we define a triangular matrix with known eigenvalues, and a random vector $b$.

```{code-cell}
ev = 10 + arange(1, 101)
A = triu(random.rand(100, 100), 1) + diag(ev)
b = random.rand(100)
```

Next we build up the first ten Krylov matrices iteratively, using renormalization after each matrix-vector multiplication.

```{code-cell}
Km = zeros([100, 30])
Km[:, 0] = b
for m in range(29):
    v = A @ Km[:, m]
    Km[:, m + 1] = v / norm(v)
```

Now we solve least-squares problems for Krylov matrices of increasing dimension, recording the residual in each case.

```{code-cell}
from numpy.linalg import lstsq
resid = zeros(30)
resid[0] = norm(b)
for m in range(1, 30):
    z = lstsq(A @ Km[:, :m], b, rcond=None)[0]
    x = Km[:, :m] @ z
    resid[m] = norm(b - A @ x)
```

The linear system approximations show smooth linear convergence at first, but the convergence stagnates after only a few digits have been found.

```{code-cell}
semilogy(range(30), resid, "-o")
xlabel("$m$"),  ylabel("$\\| b-Ax_m \\|$")
title(("Residual for linear systems"));
```
``````

(demo-subspace-arnoldi-python)=
``````{dropdown} @demo-subspace-arnoldi
:open:
We illustrate a few steps of the Arnoldi iteration for a small matrix.

```{code-cell}
A = random.choice(range(10), (6, 6))
print(A)
```

The seed vector we choose here determines the first member of the orthonormal basis.

```{code-cell}
u = random.randn(6)
Q = zeros([6, 3])
Q[:, 0] = u / norm(u)
```

Multiplication by $\mathbf{A}$ gives us a new vector in $\mathcal{K}_2$.

```{code-cell}
Aq = A @ Q[:, 0]
```

We subtract off its projection in the previous direction. The remainder is rescaled to give us the next orthonormal column.

```{code-cell}
v = Aq - dot(Q[:, 0], Aq) * Q[:, 0]
Q[:, 1] = v / norm(v)
```

On the next pass, we have to subtract off the projections in two previous directions.

```{code-cell}
Aq = A @ Q[:, 1]
v = Aq - dot(Q[:, 0], Aq) * Q[:, 0] - dot(Q[:, 1], Aq) * Q[:, 1]
Q[:, 2] = v / norm(v)
```

At every step, $\mathbf{Q}_m$ is an ONC matrix.

```{code-cell}
print(f"should be near zero: {norm(Q.T @ Q - eye(3)):.2e}")
```

And $\mathbf{Q}_m$ spans the same space as the three-dimensional Krylov matrix.

```{code-cell}
from numpy.linalg import matrix_rank
K = stack([u, A @ u, A @ A @ u], axis=-1)
Q_and_K = hstack([Q, K])
print(matrix_rank(Q_and_K))
```
``````

### 8.5 @section-krylov-gmres

(demo-gmres-intro-python)=
``````{dropdown} @demo-gmres-intro
:open:
We define a triangular matrix with known eigenvalues and a random vector $\mathbf{b}$.

```{code-cell}
ev = 10 + arange(1, 101)
A = triu(random.rand(100, 100), 1) + diag(ev)
b = random.rand(100)
```

Instead of building the Krylov matrices, we use the Arnoldi iteration to generate equivalent orthonormal vectors.

```{code-cell}
Q, H = FNC.arnoldi(A, b, 60)
print(H[:5, :5])
```

The Arnoldi bases are used to solve the least-squares problems defining the GMRES iterates.

```{code-cell}
from numpy.linalg import lstsq
resid = zeros(61)
resid[0] = norm(b)
for m in range(1, 61):
    s = hstack([norm(b), zeros(m)])
    z = lstsq(H[: m + 1, :m], s, rcond=None)[0]
    x = Q[:, :m] @ z
    resid[m] = norm(b - A @ x)
```

The approximations converge smoothly, practically all the way to machine epsilon.

```{code-cell}
semilogy(range(61), resid, "-o")
xlabel("$m$"),  ylabel("$\| b-Ax_m \|$")
title("Residual for GMRES");
```
``````

(demo-gmres-restart-python)=
``````{dropdown} @demo-gmres-restart
:open:
The following experiments are based on a matrix resulting from discretization of a partial differential equation.

```{code-cell}
d = 50;  n = d**2
A = FNC.poisson2d(d)
b = ones(n)
spy(A);
```

```{index} ! Python; gmres
```

We compare unrestarted GMRES with three different thresholds for restarting. Here we are using `gmres` from `scipy.sparse.linalg`, since our simple implementation does not offer restarting. We're also using a trick to accumulate the vector of residual norms as it runs.

```{code-cell}
from scipy.sparse.linalg import gmres
ctr = lambda rvec: resid.append(norm(rvec))
resid = [1.]
x, flag = gmres(A, b, restart=None, rtol=1e-8, atol=1e-14, maxiter=120, callback=ctr)
semilogy(resid); 
xlabel("$m$"), ylabel("residual norm")
title(("Convergence of unrestarted GMRES"));
```

```{code-cell}
maxit = 120
rtol = 1e-8
restarts = [maxit, 20, 40, 60]
hist = lambda rvec: resid.append(norm(rvec))
for r in restarts:
    resid = [1.]
    x, flag = gmres(A, b, restart=r, rtol=rtol, atol=1e-14, maxiter=maxit, callback=hist)
    semilogy(resid)

ylim(1e-8, 2)
legend(["none", "20", "40", "60"])
title(("Convergence of restarted GMRES"));
```

The "pure" GMRES curve is the lowest one. All of the other curves agree with it until the first restart. Decreasing the restart value makes the convergence per iteration generally worse, but the time required per iteration smaller as well.

``````
### 8.6 @section-krylov-minrescg

(demo-minrescg-indefinite-python)=
``````{dropdown} @demo-minrescg-indefinite
:open:
The following matrix is indefinite.

```{code-cell}
from numpy.linalg import eig
import scipy.sparse as sp
A = FNC.poisson2d(10) - 20*sp.eye(100)
ev, _ = eig(A.todense())
num_negative_ev = sum(ev < 0)
print(f"There are {sum(ev < 0)} negative and {sum(ev > 0)} positive eigenvalues")
```

We can compute the relevant quantities from {numref}`Theorem {number} <theorem-minrescg-indefinite>`.

```{code-cell}
m, M = min(-ev[ev < 0]), max(-ev[ev < 0])
kappa_minus = M / m
m, M = min(ev[ev > 0]), max(ev[ev > 0])
kappa_plus = M / m
S = sqrt(kappa_plus * kappa_minus)
rho = sqrt((S - 1) / (S + 1))
print(f"Condition numbers: {kappa_minus:.2e}, {kappa_plus:.2e}")
print(f"Convergence rate: {rho:.3f}")
```

Because the iteration number $m$ is halved in {eq}`minres-conv`, the rate constant of linear convergence is the square root of this number, which makes it even closer to 1.

Now we apply MINRES to a linear system with this matrix, and compare the observed convergence to the upper bound from the theorem.

```{index} ! Python; minres
```

```{code-cell}
from scipy.sparse.linalg import minres
b = random.rand(100)
resid = [norm(b)]
hist = lambda x: resid.append(norm(b - A @ x))
x, flag = minres(A, b, rtol=1e-8, maxiter=1000, callback=hist)
```

```{code-cell}
:tags: [hide-input]
semilogy(resid, ".-");
upper = norm(b) * rho**arange(len(resid))
semilogy(upper, "k--")
xlabel("$m$"),  ylabel("residual norm")
legend(["MINRES", "upper bound"], loc="lower left")
title("Convergence of MINRES");
```

The upper bound turns out to be pessimistic here, especially in the later iterations. However, you can see a slow linear phase in the convergence that tracks the bound closely.
``````

(demo-minrescg-converge-python)=
``````{dropdown} @demo-minrescg-converge
:open:
We will compare MINRES and CG on some quasi-random SPD problems.  The first matrix has a condition number of 100.

```{code-cell}
n = 5000
density = 0.001
A = FNC.sprandsym(n, density, rcond=1e-2)
x = arange(1, n+1) / n
b = A @ x 
```

```{index} ! Python; cg, Python; minres
```

Now we apply both methods and compare the convergence of the system residuals, using implementations imported from `scipy.sparse.linalg`.

```{code-cell}
from scipy.sparse.linalg import cg, minres
hist = lambda x: resid.append(norm(b - A @ x))

resid = [norm(b)]
xMR, flag = minres(A, b, rtol=1e-12, maxiter=100, callback=hist)
semilogy(resid / norm(b), label="MINRES")

resid = [norm(b)]
xCG, flag = cg(A, b, rtol=1e-12, maxiter=100, callback=hist)
semilogy(resid / norm(b), label="CG")

xlabel("Krylov dimension $m$"), ylabel("$\\|r_m\\| / \\|b\\|$")
grid(),  legend(),  title("Convergence of MINRES and CG");
```

There is little difference between the two methods here. Both achieve relative residual of $10^{-6}$ in aout 60 iterations, for example. The final errors are similar, too.

```{code-cell}
print(f"MINRES error: {norm(xMR - x) / norm(x):.2e}")
print(f"CG error: {norm(xCG - x) / norm(x):.2e}")
```

Next, we increase the condition number of the matrix by a factor of 25. The rule of thumb predicts that the number of iterations required should increase by a factor of about 5; i.e., 300 iterations to reach $10^{-6}$.

```{code-cell}
A = FNC.sprandsym(n, density, rcond=1e-2 / 25)
x = arange(1, n+1) / n
b = A @ x 
```

```{code-cell}
:tags: [hide-input]
from scipy.sparse.linalg import cg, minres
hist = lambda x: resid.append(norm(b - A @ x))

resid = [norm(b)]
xMR, flag = minres(A, b, rtol=1e-12, maxiter=400, callback=hist)
semilogy(resid / norm(b), label="MINRES")

resid = [norm(b)]
xCG, flag = cg(A, b, rtol=1e-12, maxiter=400, callback=hist)
semilogy(resid / norm(b), label="CG")

xlabel("Krylov dimension $m$"), ylabel("$\\|r_m\\| / \\|b\\|$")
grid(),  legend(),  title("Convergence of MINRES and CG")

print(f"MINRES error: {norm(xMR - x) / norm(x):.2e}")
print(f"CG error: {norm(xCG - x) / norm(x):.2e}")
```

Both methods have an early superlinear phase that allow them to finish slightly sooner than the factor of 5 predicted: {numref}`Theorem {number} <theorem-minrescg-converge>` is an upper bound, not necessarily an approximation. 
``````

### 8.7 @section-krylov-matrixfree

(demo-matrixfree-blur-python)=
``````{dropdown} @demo-matrixfree-blur
:open:
We use a readily available test image.

```{code-cell}
from skimage import data as testimages
from skimage.color import rgb2gray
img = getattr(testimages, "coffee")()
X = rgb2gray(img)
imshow(X, cmap="gray");
```

We define the one-dimensional tridiagonal blurring matrices.

```{code-cell}
import scipy.sparse as sp
def blurmatrix(d):
    data = [[0.25] * (d-1), [0.5] * d, [0.25] * (d-1)]
    return sp.diags(data, [-1, 0, 1], shape=(d, d))

m, n = X.shape
B = blurmatrix(m)
C = blurmatrix(n)
```

Finally, we show the results of using $k=12$ repetitions of the blur in each direction.

```{code-cell}
from scipy.sparse.linalg import matrix_power
blur = lambda X: matrix_power(B, 12) @ X @ matrix_power(C, 12)

imshow(blur(X), cmap="gray")
title("Blurred image");
```
``````

(demo-matrixfree-deblur-python)=
``````{dropdown} @demo-matrixfree-deblur
:open:
We repeat the earlier process to blur an original image $\mathbf{X}$ to get $\mathbf{Z}$.

```{code-cell}
:tags: [hide-cell]

img = getattr(testimages, "coffee")()
X = rgb2gray(img)
m, n = X.shape

import scipy.sparse as sp
def blurmatrix(d):
    data = [[0.25] * (d-1), [0.5] * d, [0.25] * (d-1)]
    return sp.diags(data, [-1, 0, 1], shape=(d, d))
B = blurmatrix(m)
C = blurmatrix(n)

from scipy.sparse.linalg import matrix_power
blur = lambda X: matrix_power(B, 12) @ X @ matrix_power(C, 12)
```

```{code-cell}
Z = blur(X)
imshow(Z, cmap="gray")
title("Blurred image");
```

Now we imagine that $\mathbf{X}$ is unknown and that we want to recover it from $\mathbf{Z}$. We first need functions that translate between vector and matrix representations.

```{code-cell}
from scipy.sparse.linalg import LinearOperator
vec = lambda Z: Z.reshape(m * n)
unvec = lambda z: z.reshape(m, n)
xform = lambda x: vec(blur(unvec(x)))
```

```{index} ! Python; LinearOperator
```

Now we declare the three-step blur transformation as a `LinearOperator`, supplying also the size of the vector form of an image.

```{code-cell}
T = LinearOperator((m * n, m * n), matvec=xform)
```

The blurring operators are symmetric, so we apply `minres` to the composite blurring transformation `T`.

```{code-cell}
from scipy.sparse.linalg import gmres
y, flag = gmres(T, vec(Z), rtol=1e-5, maxiter=50)
Y = unvec(maximum(0, minimum(1, y)))


subplot(1, 2, 1),  imshow(X, cmap="gray")
axis("off"),  title("Original")
subplot(1, 2, 2),  imshow(Y, cmap="gray")
axis("off"),  title("Deblurred");
```
``````

### 8.8 @section-krylov-precond

(demo-precond-diagonal-python)=
``````{dropdown} @demo-precond-diagonal
:open:
Here is an SPD matrix that arises from solving partial differential equations.

```{code-cell}
from scipy.sparse import sparray
import rogues
A = rogues.wathen(60, 60)
n = A.shape[0]
print(f"Matrix is {n} x {n} with {A.nnz} nonzeros")
```

```{index} ! Julia; DiagonalPreconditioner
```

There is an easy way to use the diagonal elements of $\mathbf{A}$, or any other vector, as a diagonal preconditioner.

```{code-cell}
import scipy.sparse as sp
prec = sp.diags(1 / A.diagonal(), 0)
```

We now compare CG with and without the preconditioner.

```{code-cell}
:tags: [hide-input]

from scipy.sparse.linalg import cg
b = ones(n)
hist = lambda x: resid.append(norm(b - A @ x))
resid = [norm(b)]
start = timer()
x, _ = cg(A, b, rtol=1e-4, maxiter=200, callback=hist)
print(f"No preconditioner: Finished in {timer() - start:.2f} sec")
resid_plain = resid.copy()
resid = [norm(b)]
start = timer()
x, _ = cg(A, b, rtol=1e-4, maxiter=200, M=prec, callback=hist)
print(f"Diagonal preconditioner: Finished in {timer() - start:.2f} sec")
resid_prec = resid.copy()

semilogy(resid_plain, label="no preconditioner")
semilogy(resid_prec, label="diagonal preconditioner")
xlabel("iteration"), ylabel("residual norm")
legend(),  title("Convergence of CG with and without preconditioning");
```

The diagonal preconditioner cut down substantially on the number of iterations and the execution time.
``````

(demo-precond-gmres-python)=
``````{dropdown} @demo-precond-gmres
:open:
Here is a random nonsymmetric matrix.

```{code-cell}
import scipy.sparse as sp
n = 8000
A = 2.8 * sp.eye(n) + sp.rand(n, n, 0.002)
```

Without a preconditioner, GMRES can solve a system with this matrix.

```{code-cell}
from scipy.sparse.linalg import gmres

b = random.rand(n)
hist = lambda rvec: resid.append(norm(rvec))
resid = [1.]

start = timer()
x, flag = gmres(A, b, maxiter=300, rtol=1e-10, restart=50, callback=hist)
print(f"time for plain GMRES: {timer() - start:.3f} sec")
resid_plain = resid.copy()
```

```{index} ! Python; spilu
```

The following version of incomplete LU factorization drops all sufficiently small values (i.e., replaces them with zeros). This keeps the number of nonzeros in the factors under control.

```{code-cell}
from scipy.sparse.linalg import spilu
iLU = spilu(A, drop_tol=0.2)
print(f"Factors have {iLU.nnz} nonzeros, while A has {A.nnz}")
```

The result is not a true factorization of the original matrix. However, it's close enough for an approximate inverse in a preconditioner. 

```{code-cell}
from scipy.sparse.linalg import LinearOperator
prec = LinearOperator((n, n), matvec=lambda y: iLU.solve(y))

resid = [1.];  start = timer()
x, flag = gmres(A, b, M=prec, maxiter=300, rtol=1e-10, restart=50, callback=hist)
print(f"time for preconditioned GMRES: {timer() - start:.3f} sec")
resid_prec = resid
```

```{code-cell}
:tags: [hide-input]
semilogy(resid_plain, label="no prec.")
semilogy(resid_prec, label="iLU prec.")
xlabel("iteration number"),  ylabel("residual norm")
legend()
title("GMRES convergence compared");
``` 
``````