---
numbering:
  enumerator: 8.1.%s
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

(section-krylov-structure)=

# Sparsity and structure

```{index} ! sparse matrix
```

Very large matrices cannot be stored all within primary memory of a computer unless they are **sparse**. A sparse matrix has *structural zeros*, meaning entries that are known to be exactly zero. For instance, the adjacency matrix of a graph has zeros where there are no links in the graph. To store and operate with a sparse matrix efficiently, it is not represented as an array of all of its values. There is a variety of sparse formats available; for the most part, you can imagine that the matrix is stored as triples $(i,j,A_{ij})$ for all the nonzero $(i,j)$ locations.

## Computing with sparse matrices

```{index} adjacency matrix
```

Most graphs with real applications have many fewer edges than the maximum possible $n^2$ for $n$ nodes. Accordingly, their adjacency matrices have mostly zero elements and should be represented sparsely.

::::{prf:example} Sparsity
:label: demo-structure-sparse

```{tip}
Functions to work with sparse matrices are found in the `scipy.sparse` module.
```

Here we load the adjacency matrix of a graph with 2790 nodes from a [data file](https://raw.github.com/fncbook/fnc/master/python/roswelladj.mat). Each node is a web page referring to Roswell, NM, and the edges represent links between web pages. (Credit goes to Panayiotis Tsaparas and the University of Toronto for making this data public.)

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

::::

```{index} fill-in of sparse matrices
```

Arithmetic operations such as `+`, `-`, `*`, and `^` respect and exploit sparsity if the matrix operands are sparse. However, matrix operations may substantially decrease the amount of sparsity, a phenomenon known as **fill-in**.

::::{prf:example} Fill-in of a sparse matrix
:label: demo-structure-fill


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

::::

## Banded matrices

```{index} banded matrix
```

A particularly important type of sparse matrix is a banded matrix. Recall from {numref}`section-linsys-structure` that $\mathbf{A}$ has **upper bandwidth** $p$ if $j-i > p$ implies $A_{ij}=0$, and **lower bandwidth** $q$ if $i-j > q$ implies $A_{ij}=0$. We say the total **bandwidth** is $p+q+1$. Banded matrices appear naturally in many applications where each element interacts directly with only a few neighbors.

```{index} LU factorization
```

Without pivoting, an LU factorization preserves bandwidth, but pivoting can change or destroy bandedness.

::::{prf:example} Banded matrices
:label: demo-structure-sparseband

```{index} ! Python; diags
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

::::

## Systems and eigenvalues

If given a sparse matrix, the backslash operator will automatically try a form of sparse-aware Cholesky or pivoted LU factorization. Depending on the sparsity pattern of the matrix, the time taken to solve the linear system may be well below the $O(n^3)$ needed in the general case.

```{index} eigenvalue decomposition
```

For very large matrices, it's unlikely that you will want to find all of its eigenvalues and eigenvectors. In {numref}`section-krylov-subspace` we describe some of the math behind an algorithm that can find a selected number of eigenvalues of largest magnitude, lying to the extreme left or right, or nearest a given complex number.

::::{prf:example} Eigenvalues of sparse matrices
:label: demo-structure-linalg


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
from scipy.linalg import solve
F = A.todense()
start = timer()
xx = solve(F, b)
print(f"dense time: {timer() - start:.3g} sec")
print(f"residual: {norm(b - A @ xx, 2):.1e}")
```

::::

## Exercises

``````{exercise}
:label: problem-structure-spdiagm
⌨ Use `spdiagm` to build the $50\times 50$ matrices

```{math}
\mathbf{A} =
\begin{bmatrix}
-2 & 1 & & &  \\
1 & -2 & 1 & &  \\
& \ddots & \ddots & \ddots & \\
& & 1 & -2 & 1 \\
& & & 1 & -2
\end{bmatrix}, \qquad
\mathbf{B} =
\begin{bmatrix}
-2 & 1 & & & 1 \\
1 & -2 & 1 & &  \\
& \ddots & \ddots & \ddots & \\
& & 1 & -2 & 1 \\
1 & & & 1 & -2
\end{bmatrix}.
```

For each matrix, use `spy` and an inspection of the $5\times 5$ submatrices in the corners to verify the correctness of your matrices.
``````

``````{exercise}
:label: problem-structure-factoring
⌨ This problem requires the matrix used in @demo-structure-fill.
Download [a data file](https://raw.github.com/fncbook/fnc/master/chapter8/smallworld.mtx) by right-clicking the link and selecting "Save As...".
``` python
import scipy.io as spio
A = spio.mmread("smallworld.mtx")
```

**(a)** Find the density of $\mathbf{A}$ (number of nonzeros divided by total number of elements), $\mathbf{A}^2$, $\mathbf{A}^4$, and $\mathbf{A}^8$. (You should find that it increases with the power of $\mathbf{A}$.)

**(b)** The LU factors tend to at least partially retain sparsity. Find the density of the $\mathbf{L}$ and $\mathbf{U}$ factors of $\mathbf{A}$ using `lufact` (@function-lufact). (If you get an error, convert the matrix to dense form first.)

**(c)** Repeat part (b) for the QR factorization using `qrfact` (@function-qrfact). (If you get an error, convert the matrix to dense form first.)
``````

``````{exercise}
:label: problem-structure-roswell
⌨ One use of adjacency matrices is to analyze the links between members of a collection. Obtain the adjacency matrix $\mathbf{A}$ from @demo-structure-sparse via the following:

Download [roswell.mtx](roswell.mtx) by clicking the link and saving (you may need to fix the file name).
``` python
import scipy.io as spio
A = spio.mmread("roswell.mtx")
```

The matrix catalogs the links between web sites related to the town of Roswell, NM, with $A_{ij}=1$ if and only if site $i$ links to site $j$.

**(a)** Verify numerically that the matrix does not include any links from a site to itself.

**(b)** Verify numerically that $\mathbf{A}$ is not symmetric. (Thus, its graph is a directed one.)

**(c)** How many sites in the group are not pointed to by any other sites in the group?

**(d)** Which site points to the most other sites?

**(e)** Which site is pointed to the most by the other sites? This is a crude way to establish the most important site.

**(f)** There are $2790^2$ possible ways to connect ordered pairs of sites. What fraction of these pairs is connected by a walk of links that is no greater than three in length?


```{index} ! graph Laplacian matrix, ! degree matrix
```
``````

``````{exercise}
:label: problem-structure-laplacian
⌨ The [graph Laplacian matrix](wiki:Laplacian_matrix) is $\mathbf{L}=\mathbf{D}-\mathbf{A}$, where $\mathbf{A}$ is the adjacency matrix and $\mathbf{D}$ is the *degree matrix*, a diagonal matrix with diagonal entries $d_{jj}=\sum_{i=1}^n a_{ij}$. 

Follow the directions in @problem-structure-roswell to obtain an adjacency matrix $\mathbf{A}$. Then find the five eigenvalues of $\mathbf{L}$ having largest magnitude.
``````

``````{exercise}
:label: problem-structure-actorsmat
⌨ See @problem-insight-actors for instructions on loading a matrix $\mathbf{A}$ that contains information about the appearances of 392,400 actors in 127,823 movies, as given by the [Internet Movie Database](wiki:IMDb). Specifically, $A_{ij}=1$ if actor $j$ appeared in movie $i$, and all other elements are zero.

**(a)** What is the maximum number of actors appearing in any one movie?

**(b)** How many actors appeared in exactly three movies?

**(c)** Define $\mathbf{C}=\mathbf{A}^T\mathbf{A}$. How many nonzero entries does $\mathbf{C}$ have? What is the interpretation of $C_{ij}$?
``````

```{index} ! Helmholtz equation
```

``````{exercise}
:label: problem-structure-helmholtz
⌨  A matrix that arises from the [Helmholtz equation](wiki:Helmholtz_equation) for wave propagation can be specified using 

```julia
A = FNC.poisson(n) - k^2*I;
```
where $k$ is a real parameter. Let $n=50$. 

**(a)** Let $k=1$. What is the size of $\mathbf{A}$? What is its density?

**(b)** Still with $k=1$, use `eigs` to find the four largest and four smallest (in magnitude) eigenvalues of $\mathbf{A}$. (See @demo-structure-linalg for examples.)

**(c)** The eigenvalues are all real. Find a value of $k$ so that $\mathbf{A}$ has exactly three negative eigenvalues.
``````
