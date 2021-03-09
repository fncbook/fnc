# Sparsity and structure

```{index} matrix; sparse
``` 
```{index} sparse matrix
```

Very large matrices cannot be stored all within primary memory of a computer unless they are {term}`sparse`. A sparse matrix has *structural zeros*, meaning entries that are known to be exactly zero.} For instance, the adjacency matrix of a graph has zeros where there are no links in the graph. To store and operate with a sparse matrix efficiently, it is not represented as an array of all of its values. There is a variety of sparse formats available; MATLAB chooses *compressed sparse column* format. For the most part, you can imagine that the matrix is stored as triples $(i,j,A_{ij})$ for all the nonzero $(i,j)$ locations.

## Computing with sparse matrices

```{index} adjacency matrix
```

Most graphs with real applications have many fewer edges than the maximum possible $n^2$ for $n$ nodes. Accordingly, their adjacency matrices have mostly zero elements and should be represented sparsely. 

::::{prf:example} Julia demo
:class: demo
:label: demos-structure-sparse
{doc}`demos/structure-sparse`
::::

Arithmetic operations such as `+`, `-`, `*`, and `^` respect and exploit sparsity, if the matrix operands are sparse. However, matrix operations may substantially decrease the amount of sparsity, a phenomenon known as {term}`fill-in`.

In the case of an adjacency matrix $\mathbf{A}$, for example, the $(i,j)$ entry of matrix $\mathbf{A}^k$ for positive integer $k$ is the number of paths of length $k$ from node $i$ to node $j$.


::::{prf:example} Julia demo
:class: demo
:label: demos-structure-fill
{doc}`demos/structure-fill`
::::

## Banded matrices

```{index} matrix; banded
```

A particularly important type of sparse matrix is a banded matrix. [Earlier](../linsys/structure.md) we said that $\mathbf{A}$ has **upper bandwith** $p$ if $j-i > p$ implies $A_{ij}=0$, and **lower bandwidth** $q$ if $i-j > q$ implies $A_{ij}=0$. We say the total **bandwidth** is $p+q+1$. 

Banded matrices appear naturally in many applications where each node interacts directly with only its closest neighbors. Without pivoting, an LU factorization preserves bandwidth, but pivoting can change or destroy bandedness.

```{index} matrix; factorization
```

::::{prf:example} Julia demo
:class: demo
:label: demos-structure-banded
{doc}`demos/structure-banded`
::::

## Linear systems and eigenvalues

If given a sparse matrix, the backslash operator will automatically try a form of sparse-aware Cholesky or pivoted LU factorization. Depending on the sparsity pattern of the matrix, the time taken to solve the linear system may be well below the $O(n^3)$ needed in the general case.

```{index} eigenvalue decomposition
```

For very large matrices, it's unlikely that you will want to find all of its eigenvalues and eigenvectors. In an [upcoming section](subspace.md) we will describe some of the math behind an algorithm that can find a selected number of eigenvalues of largest magnitude, lying to the extreme left or right, or nearest a given complex number. 

::::{prf:example} Julia demo
:class: demo
:label: demos-structure-linalg
{doc}`demos/structure-linalg`
::::

1. ⌨ Use `spdiagm` to build the $200\times 200$ matrices
  
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

    ::::{only} solutions
    %%
    d = repmat([1 -2 1],200,1);
    A = spdiags(d,-1:1,200,200);
    spy(A)
    full( A(1:5,1:5) )

    %%
    % This matrix is the same, except for corner elements.
    d = repmat([1 1 -2 1 1],200,1);
    B = spdiags(d,[-199 -1 0 1 199],200,200);
    spy(B)
    full( B(1:5,end-4:end) )
    ::::

2. ⌨ This problem uses the `smallworld.jld2` matrix as seen in {prf:ref}`demos-structure-fill`. 
  
    **(a)** Find the density of $\mathbf{A}$ (number of nonzeros divided by total number of elements), $\mathbf{A}^2$, $\mathbf{A}^4$, and $\mathbf{A}^8$. (You should find that it increases with the power of $A$.)

    **(b)** The LU factors tend to at least partially retain sparsity. Find the density of the $\mathbf{L}$ and $\mathbf{U}$ factors of $\mathbf{A}$ using the built-in `lu` (with two outputs).

    **(c)** The QR factors have no particular relationship to sparsity in general. Repeat part (b) for the QR factorization using `qr`.

    ::::{only} solutions
    %%
    A = bucky;
    N = numel(A);
    %%
    den1 = nnz(A)/N
    den2 = nnz(A^2)/N
    den4 = nnz(A^4)/N
    den8 = nnz(A^8)/N
    %%
    [L,U] = lu(A);
    denL = nnz(L)/N
    denU = nnz(U)/N
    %%
    [Q,R] = qr(A);
    denQ = nnz(Q)/N
    denR = nnz(R)/N
    ::::

3. not available

4. ✍ Prove the statement illustrated in {prf:ref}`demos-structure-fill`: if $\mathbf{A}$ is an adjacency matrix, the $(i,j)$ entry of matrix $\mathbf{A}^k$ for positive integer $k$ is the number of paths of length $k$ from node $i$ to node $j$.

5. ⌨ One use of adjacency matrices is to analyze the links between members of a collection. Obtain an adjacency matrix $\mathbf{A}$ by loading the file `roswelladj.jld2` as in {prf:ref}`demos-structure-sparse`. (Credit goes to Panayiotis Tsaparas and the University of Toronto for making this data public.) The matrix catalogs the links between web sites related to the town of Roswell, NM, with $A_{ij}=1$ if and only if site $i$ points to site $j$.
  
    **(a)** Verify numerically that the matrix does not include any links from a site to itself.

    **(b)** Note that $\mathbf{A}$ does not need to be symmetric, because web links point in one direction. Verify numerically that $\mathbf{A}$ is not symmetric.

    **(c)** How many sites in the group are not pointed to by any other sites in the group?

    **(d)** Which site points to the most other sites?

    **(e)** Which site is pointed to the most by the other sites? This is a crude way to establish the most important site.

    **(f)** There are $2790^2$ possible ways to connect ordered pairs of sites. What fraction of these pairs is connected by a path of links that is no greater than three in length?

    ::::{only} solutions
    %%
    %%
    has_self_links = norm( diag(A), inf )   % no self-links

    %%
    S = A - A';
    is_symmetric = norm( S(:), inf )      % not symmetric

    %%
    % How many are not pointed to?
    not_pointed_to = find( sum(A,1)==0 );
    num_not_pointed_to = length(not_pointed_to)

    %%
    % Which site points to the most others?
    num_pointed_to = sum(A,2);
    [~,points_to_most] = max(num_pointed_to)

    %%
    % Which site is most pointed at?
    num_pointed_at = sum(A,1);
    [~,pointed_at_most] = max(num_pointed_at)

    %%
    % What fraction of pairs are connected by paths of length <=3?
    A2 = A^2;   % # paths of length 2
    A3 = A^3;   % # paths of length 3
    S = A + A2 + A3;
    fraction_of_paths = nnz(S) / 2790^2
    ::::

    ```{index} matrix; graph Laplacian
    ```
    ```{index} matrix; degree
    ```

6. ⌨ The *graph Laplacian matrix* is $\mathbf{L}=\mathbf{D}-\mathbf{A}$, where $\mathbf{A}$ is the adjacency matrix and $\mathbf{D}$ is the *degree matrix*, a diagonal matrix with diagonal entries $d_{jj}=\sum_{i=1}^n a_{ij}$. Load the file `roswelladj.jld2`, available from the book website, to get an adjacency matrix $\mathbf{A}$. Then find the five eigenvalues of $\mathbf{L}$ having largest magnitude.

    ::::{only} solutions
    load roswelladj
    D = spdiags(sum(A)',0,2790,2790);
    L = D - A;
    eigs(L,5,'lr')
    ::::

7. ⌨ \label{pro:actorsmat} The file `actors.jld2` is available at the book site.  Based on data provided by the Self-Organized Networks Database at the University of Notre Dame, it contains information about the appearances of 392,400 actors in 127,823 movies, as given by the Internet Movie Database. The matrix $\mathbf{A}$ has $A_{ij}=1$ if actor $j$ appeared in movie $i$ and zero elsewhere.

    **(a)** What is the maximum number of actors appearing in any one movie?

    **(b)** How many actors appeared in exactly three movies?

    **(c)** Define $\mathbf{C}=\mathbf{A}^T\mathbf{A}$. How many nonzero entries does $\mathbf{C}$ have? What is the interpretation of $C_{ij}$?


