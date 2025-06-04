---
numbering:
  enumerator: 8.1.%s
---
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
`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-structure-sparse-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-structure-sparse-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-structure-sparse-python
:::
````
`````
::::

```{index} fill-in of sparse matrices
```

Arithmetic operations such as `+`, `-`, `*`, and `^` respect and exploit sparsity if the matrix operands are sparse. However, matrix operations may substantially decrease the amount of sparsity, a phenomenon known as **fill-in**.

::::{prf:example} Fill-in of a sparse matrix
:label: demo-structure-fill
`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-structure-fill-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-structure-fill-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-structure-fill-python
:::
````
`````
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
`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-structure-sparseband-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-structure-sparseband-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-structure-sparseband-python
:::
````
`````
::::


## Systems and eigenvalues

If given a sparse matrix, the backslash operator will automatically try a form of sparse-aware Cholesky or pivoted LU factorization. Depending on the sparsity pattern of the matrix, the time taken to solve the linear system may be well below the $O(n^3)$ needed in the general case.

```{index} eigenvalue decomposition
```

For very large matrices, it's unlikely that you will want to find all of its eigenvalues and eigenvectors. In {numref}`section-krylov-subspace` we describe some of the math behind an algorithm that can find a selected number of eigenvalues of largest magnitude, lying to the extreme left or right, or nearest a given complex number. 

::::{prf:example} Eigenvalues of sparse matrices
:label: demo-structure-linalg
`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-structure-linalg-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-structure-linalg-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-structure-linalg-python
:::
````
`````
::::

## Exercises

``````{exercise}
:label: problem-structure-spdiagm
⌨ Use `spdiagm` to build the $50\times 50$ matrices

:::{math}
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
:::

For each matrix, use `spy` and an inspection of the $5\times 5$ submatrices in the corners to verify the correctness of your matrices.
``````

``````{exercise}
:label: problem-structure-factoring
⌨ This problem requires the matrix used in @demo-structure-fill.
`````{tab-set} 
````{tab-item} Julia
:sync: julia
Download [smallworld.mat](smallworld.mat) by clicking the link and saving (you may need to fix the file name).
``` julia
using MAT
A = matread("smallworld.mat")["A"]
```
```` 

````{tab-item} MATLAB
:sync: matlab
Download [smallworld.mat](smallworld.mat) by clicking the link and saving (you may need to fix the file name).
``` matlab
load smallworld
```
```` 

````{tab-item} Python
:sync: python
Download [smallworld.mtx](smallworld.mtx) by clicking the link and saving (you may need to fix the file name).
``` python
import scipy.io as spio
A = spio.mmread("smallworld.mtx")
```
```` 
`````

**(a)** Find the density of $\mathbf{A}$ (number of nonzeros divided by total number of elements), $\mathbf{A}^2$, $\mathbf{A}^4$, and $\mathbf{A}^8$. (You should find that it increases with the power of $\mathbf{A}$.)

**(b)** The LU factors tend to at least partially retain sparsity. Find the density of the $\mathbf{L}$ and $\mathbf{U}$ factors of $\mathbf{A}$ using `lufact` (@function-lufact). (If you get an error, convert the matrix to dense form first.)

**(c)** Repeat part (b) for the QR factorization using `qrfact` (@function-qrfact). (If you get an error, convert the matrix to dense form first.)
``````

``````{exercise}
:label: problem-structure-roswell
⌨ One use of adjacency matrices is to analyze the links between members of a collection. Obtain the adjacency matrix $\mathbf{A}$ from @demo-structure-sparse via the following:

`````{tab-set} 
````{tab-item} Julia
:sync: julia
Download [roswell.mat](roswell.mat) by clicking the link and saving (you may need to fix the file name).
``` julia
using MAT
A = matread("roswell.mat")["A"]
```
```` 

````{tab-item} MATLAB
:sync: matlab
Download [roswell.mat](roswell.mat) by clicking the link and saving (you may need to fix the file name).
``` matlab
load roswell
```
```` 

````{tab-item} Python
:sync: python
Download [roswell.mtx](roswell.mtx) by clicking the link and saving (you may need to fix the file name).
``` python
import scipy.io as spio
A = spio.mmread("roswell.mtx")
```
```` 
`````

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
