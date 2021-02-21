\chapter{Matrix analysis}
\label{cha:matrixanaly}

\graphicspath{{matrixanaly/figures/}{matrixanaly/examples/}}

\begin{flushright}
  \small \sf
  \begin{tabular}{p{3in}}

    Judge me by my size, do you?\\ \hfill --- Yoda, *The Empire Strikes Back* 
```{index} Yoda
```
 
```{index} The Empire Strikes Back
```

  \end{tabular}
\end{flushright}


In previous chapters we have seen how matrices that represent square or overdetermined linear systems of equations can be manipulated into LU and QR factorizations. But matrices have other factorizations that are more intrinsic to their nature as mathematical linear transformations. The most fundamental of these are the eigenvalue and singular value decompositions.

These decompositions can be used to solve linear and least squares systems, but they have greater value in how they represent the matrix itself. They lead to critical and quantitative insights about the structure of the underlying transformation and suggest ways to approximate it efficiently. In this chapter we will look at both of these fundamental decompositions and hint at just a few of their computational applications.

\newpage

# From matrix to insight
\label{sec:matrix-insight}

\texthighlight{matrixinsight}{Any two-dimensional array of numbers may be interpreted as a matrix.} Whether or not this is the only point of view that matters to a particular application, it does lead to particular types of analysis. The related mathematical and computational tools are universally applicable and find diverse uses.

## Tables as matrices

Tables are used to represent variation of a quantity with respect to two variables. These variables may be encoded as the rows and columns of a matrix.


::::{proof:example}
  Suppose we have a *corpus*, or collection of text documents. A *term-document matrix* has one column for each document and one row for each unique term appearing in the corpus. The $(i,j)$ entry of the matrix is the number of times term $i$ appears in document $j$. That is, column $j$ of the matrix is a term-frequency vector quantifying all occurrences of the indexed terms. A new document could be represented by its term-frequency vector, which is then comparable to the columns of the matrix. Or, a new term could be represented by counting its appearances in all of the documents and be compared to the rows of the matrix.

  It turns out that by finding the singular value decomposition (\secref{svd}) of the term-document matrix, the strongest patterns within the corpus can be isolated, frequently corresponding to what we interpret as textual meaning. This is known as *latent semantic analysis.*
::::



::::{proof:example}
  The website \url{www.congress.gov/roll-call-votes} offers data on all the votes cast in each session of the U.S.\ Congress. We can put members of Congress along the columns of a matrix and bills along the rows, recording a number that codes for "yea,"  "nay," "none," etc. The singular value decomposition (\secref{svd}) can reveal an objective, reproducible analysis of the partisanship and cooperation of individual members.
::::



::::{proof:example}
  In 2006 the online video service Netflix started an open competition for a \$1 million prize. They provided a data set of 100,480,507 ratings (one to five stars) made by 480,189 users for 17,770 movies. Each rating is implicitly an entry in a 17,770$\times$480,189 matrix. The object of the prize was to predict a user's ratings for movies they had not rated. This is known as a *matrix completion problem.* (It took 6 days for a contestant to improve on Netflix's private algorithm, and in 2009 the million-dollar prize was awarded to a team that had improved the performance by over 10\%.)
::::


## Graphs as matrices

An important concept in modern mathematics is that of a **graph**. A graph consists of a set $V$ of 
```{index} graph nodes and edges
```
 **nodes** and a set $E$ of **edges**, each of which is an ordered pair of nodes. The natural interpretation is that the edge $(v_i,v_j)$ denotes a link from node $i$ to node $j$, in which case we say that node $i$ is adjacent to node $j$. One usually visualizes small graphs by drawing points for nodes and arrows or lines for the edges.

Graphs are useful because they are the simplest way to represent link structure---of social networks, airline routes, power grids, sports teams, and web pages, to name a few examples. They also have close ties to linear algebra. Prominent among these is the 
```{index} adjacency matrix
```
 
```{index} matrix!adjacency|see {adjacency matrix
```
} {term}`adjacency matrix` of the graph. If the graph has $n$ nodes, then its $n\times n$ adjacency matrix $\mathbf{A}$ has elements

```{math}
:label: adjmat
A_{ij} =
\begin{cases}
1,& \text{if $(v_i,v_j)\in E$ (i.e., there is an edge from node $i$ to node $j$)},\\
0,& \text{otherwise}.
\end{cases}
```



::::{proof:example}
    \inputexample{matrixanaly}{buckyplot}
::::



The representation of a graph by its adjacency matrix opens up the possibility for many kinds of analysis of the graph. One might ask whether the nodes admit a natural partition into clusters, for example. Or one might ask to rank the nodes in order of importance to the network as determined by some objective criteria---an application made famous by Google's PageRank algorithm, and one which is mathematically stated as an eigenvalue problem (\secref{evd}).


## Images as matrices
%\label{sec:images}


```{index} matrix!as image
```

Computers most often represent images as rectangular arrays of pixels, each of which is colored according to numerical values for red (R), green (G), and blue (B) components of white light. Typically these are given as integers in the range from zero (no color) to 255 (full color). Thus, an image that is $m$-by-$n$ pixels can be stored as an $m$-by-$n$-by-3 array of integer values.

We will simplify the representation by considering images represented using pixels that can take only shades of gray. We will also use floating-point numbers rather than integers, so that we can operate on them using real arithmetic, though we will stay with the convention that values should be in the range $[0,255]$. (Pixels below zero or above 255 will be colored pure black or pure white, respectively.)


::::{proof:example}
    \inputexample{matrixanaly}{imreaddemo}
::::


Representation of an image as a matrix allows us to describe some common image operations in terms of linear algebra. Furthermore, the singular value decomposition (\secref{svd}) can be used to compress the information.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\begin{exercises}
  \input{matrixanaly/exercises/Insight}
\end{exercises}

\clearpage

# Eigenvalue decomposition
\label{sec:evd}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

To this point we have dealt frequently with the solution of the linear system $\mathbf{A} \mathbf{x}=\mathbf{b}$. Alongside this problem in its importance to linear algebra is the eigenvalue problem,

```{math}
  :label: eigdef
  \mathbf{A} \mathbf{x} = \lambda \mathbf{x},
```

for a scalar 
```{index} eigenvalue
```
 **eigenvalue** $\lambda$ and an associated nonzero **eigenvector** $\mathbf{x}$.

## Complex matrices

A matrix with real entries can have complex eigenvalues. Therefore we assume all matrices, vectors, and scalars may be complex in what follows. Recall that a complex number can be represented as $a+ib$, for real $a$ and $b$ and where $i^2=-1$. The **complex conjugate** of $x=a+ib$ is denoted $\bar{x}$ and is given by $\bar{x}=a-ib$. The {term}`hermitian` or conjugate transpose of a matrix $\mathbf{A}$ is denoted $\mathbf{A}^*$ and is given by $\mathbf{A}^*=(\overline{\mathbf{A}})^T=\overline{\mathbf{A}^T}$.

For the most part, "hermitian" replaces "transpose" when dealing with complex matrices. For instance, the inner product of complex vectors $\mathbf{u}$ and $\mathbf{v}$ is

```{math}
  :label: complexinnerprod
  \mathbf{u}^* \mathbf{v} = \sum_{k=1}^n \overline{u}_k v_k,
```

which in turn defines the 2-norm for complex vectors, and thereby
matrices as well. The definitions of orthogonal and orthonormal sets
of complex-valued vectors use ${}^*$
instead of ${}^T$. The analog of an orthgonal matrix in the complex
case---that is, a square matrix whose columns are orthonormal in the
complex sense---is said to be 
```{index} unitary matrix
```
 
```{index} matrix!unitary|see {unitary matrix
```
} 
```{index} orthogonal!matrix|seealso {unitary matrix
```
} {term}`unitary`. A unitary matrix $\mathbf{U}$
satisfies $\mathbf{U}^{-1}=\mathbf{U}^*$ and $\| \mathbf{U}\mathbf{x} \|_2=\| \mathbf{x} \|_2$ for any
complex vector $\mathbf{x}\in\mathbb{C}^n$.

% ## Independence and basis

% Recall that a linear combination of vectors $\mathbf{v}_1$, \ldots, $\mathbf{v}_k$
% is a vector of the form
% 
```{math}
%   c_1 \mathbf{v}_1 +  \dots + c_k \mathbf{v}_k,
% ```

% for scalars $c_1,\dots,c_k$. Linear combinations are equivalent to
% matrix-vector products, via
% The **span** of these vectors is the set of all possible linear
% combinations; i.~e., the range of $\mathbf{V}$. The vectors are said to be
% **independent** if
% 
```{math}
%   c_1 \mathbf{v}_1+\dots + c_k \mathbf{v}_k = oldsymbol{0} \text{ implies }
%   c_1=c_2 = \dots = c_k = 0.
% ```

% Otherwise they are said to be **dependent**. Dependence implies
% that the equation $V\mathbf{x}=oldsymbol{0}$ has a nonzero solution for
% $\mathbf{x}$. If $\mathbf{V}$ is square, this fact is equivalent to the singularity
% of $\mathbf{V}$.

% A **basis** of $\mathbb{C}^n$ is a set of $n$ linearly independent
% vectors. If $\{\mathbf{v}_{1},\dots,\mathbf{v}_{n}\}$ is such a set, then
% $V\mathbf{x}=\mathbf{b}$ has a solution for every $\mathbf{b}$. If $\mathbf{V}$ is square, then the
% statements "$\mathbf{V}$ is nonsingular" and ``the columns of $\mathbf{V}$ form a
% basis'' are equivalent.

% %A basis is like a coordinate system for $\mathbb{C}^n$.  If $\mathbf{V}$ is
% %square and nonsingular, its columns form a basis and the vector $y=Vx$
% %is the linear combination whose coefficients are the components of
% %$x$. In other words, the entries of $x=V^{-1}y$ are the coefficients
% %of $y$ when it is expressed in the basis made up of the columns of
% %$\mathbf{V}$. This connection between matrix-vector multiplication and changes
% %of basis is quite important.


## Eigenvalue decomposition

The eigenvalue equation $\mathbf{A}\mathbf{x}=\lambda\mathbf{x}$ is equivalent to $(\lambda\mathbf{I} - \mathbf{A})\mathbf{x}=oldsymbol{0}$, which, since $\mathbf{x}$ is nonzero in order to be an eigenvector, implies that $\lambda\mathbf{I} - \mathbf{A}$ is a singular matrix. This observation leads to the familiar property of an eigenvalue being a root of the 
```{index} characteristic polynomial
```
 **characteristic polynomial** $\det(\lambda \mathbf{I} - \mathbf{A})$. From here one concludes that an $n\times n$ matrix has $n$ eigenvalues, counting multiplicity.

Hence suppose that $\mathbf{A}\mathbf{v}_k=\lambda_k\mathbf{v}_k$ for $k=1\ldots,n$. We can summarize these as
\begin{align}
  \begin{bmatrix}
    \mathbf{A}\mathbf{v}_1 & \mathbf{A}\mathbf{v}_2 & \cdots & \mathbf{A}\mathbf{v}_n
  \end{bmatrix}
  &=
    \begin{bmatrix}
      \lambda_1 \mathbf{v}_1 & \lambda_2\mathbf{v}_2 & \cdots & \lambda_n \mathbf{v}_n
    \end{bmatrix}\notag\\
  \mathbf{A} \begin{bmatrix}
    \mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n
  \end{bmatrix}
  &=
\begin{bmatrix}
    \mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n
  \end{bmatrix}
  \begin{bmatrix}
    \lambda_1 & &  &  \\
    & \lambda_2 & & \\
    & & \ddots & \\
    & & & \lambda_n
  \end{bmatrix} \notag\\
  \mathbf{A} \mathbf{V} &= \mathbf{V} \mathbf{D}. :label: ev-all
\end{align}
So far $\mathbf{A}$ could be any square matrix. But if we also assume that $\mathbf{V}$ is a nonsingular matrix, we can rewrite {eq}`ev-all` as

```{math}
  :label: evdecomp
  \mathbf{A} = \mathbf{V} \mathbf{D} \mathbf{V}^{-1}.
```

\texthighlight{evd}{Equation {eq}`evdecomp` is called an 
```{index} matrix!factorization!|seealso {eigenvalue decomposition
```
} 
```{index} eigenvalue decomposition
```
 {term}`eigenvalue decomposition` (EVD) of $\mathbf{A}$.} If $\mathbf{A}$ has an EVD, we say that $\mathbf{A}$ is 
```{index} matrix!diagonalizable
```
 {term}`diagonalizable`; otherwise $\mathbf{A}$ is **nondiagonalizable** (or *defective*). One simple example of a nondiagonalizable matrix is

```{math}
  :label: jordanblock
  \mathbf{B} = \begin{bmatrix}
    1 & 1\\0 & 1
  \end{bmatrix}.
```

There is a common circumstance in which we can guarantee an EVD exists; the proof of the following theorem can be found in many elementary texts on linear algebra.

::::proof:theorem}
  If the $n\times n$ matrix $\mathbf{A}$ has $n$ distinct eigenvalues, then $\mathbf{A}$
  is diagonalizable.
::::


The {term}`eig` command can be used to compute the eigenvalue decomposition of a given matrix.


::::{proof:example}
  \inputexample{matrixanaly}{eigintro}
::::


Observe that if $\mathbf{A}\mathbf{v} = \lambda \mathbf{v}$ for nonzero $\mathbf{v}$, then the equation remains true for any nonzero multiple of $\mathbf{v}$, so eigenvectors are not unique, and neither is an EVD.



## Similarity and change of basis

The particular relationship between matrices $\mathbf{A}$ and $\mathbf{D}$ in {eq}`evdecomp` is important.  If $\mathbf{S}$ is any nonsingular matrix, we say that $\mathbf{B}=\mathbf{S}\mA\mathbf{S}^{-1}$ is 
```{index} matrix!similar
```
 **similar** to $\mathbf{A}$. A similarity transformation does not change eigenvalues, a fact that is typically proved in elementary linear algebra texts.
%A similarity transformation is really a change of
%basis. Interpret the product $SAS^{-1}x$ from right to left as follows:
%First express $x$ in the basis of columns of $\mathbf{S}$, then apply $\mathbf{A}$, then
%undo the change of basis.
%The following theorem states that eigenvalues are independent of basis.

::::proof:theorem}
  If $\mathbf{X}$ is an nonsingular matrix, then $\mathbf{X}\mA\mathbf{X}^{-1}$ has the same
  eigenvalues as $\mathbf{A}$.
::::


Similarity transformation has a relatively simple interpretation. First, consider the product of a nonsingular $\mathbf{X}$ with any vector:

```{math}
  \mathbf{y} = \mathbf{X} \mathbf{z} = z_1 \mathbf{x}_1 +  \dots + z_n \mathbf{x}_n.
```

We call $z_1,\ldots,z_n$ the *coordinates* of the vector $\mathbf{y}$ with respect to the columns of $\mathbf{X}$. That is, $\mathbf{z}$ is a representation of $\mathbf{y}$ relative to the basis implied by the columns of $\mathbf{X}$. But then $\mathbf{z} = \mathbf{X}^{-1} \mathbf{y}$, so left-multiplication by $\mathbf{X}^{-1}$ converts the vector $\mathbf{y}$ into that representation. In other words, \texthighlight{changebasis}{multiplication by the inverse of a matrix performs a **change of basis** into the coordinates associated with the matrix.}

So now consider the EVD {eq}`evdecomp` and the product $\mathbf{u} =\mathbf{A} \mathbf{x}$, or $(\mathbf{V}^{-1}\mathbf{u}) = \mathbf{D}(\mathbf{V}^{-1}\mathbf{x})$. This equation says that if you express the input $\mathbf{x}$ and the output $\mathbf{u}$ into the coordinates of the $\mathbf{V}$-basis, then the relationship between them
is diagonal. That is, the EVD is about finding a basis for $\mathbb{C}^n$ in which the map $\mathbf{x}\mapsto\mathbf{A}\mathbf{x}$ is a diagonal one in which the coordinates are independently rescaled.

The fact that the EVD represents a change of basis in both the domain and range spaces makes it useful for matrix powers: $$\mathbf{A}^2=(\mathbf{V}\mD\mathbf{V}^{-1})(\mathbf{V}\mD\mathbf{V}^{-1})=\mathbf{V}\mD(\mathbf{V}^{-1}\mathbf{V})\mathbf{D}\mV^{-1}=\mathbf{V}\mD^2\mathbf{V}^{-1},$$ and so on. Because $\mathbf{D}$ is diagonal, its power $\mathbf{D}^k$ is just the diagonal matrix of the powers of the eigenvalues.


## Conditioning of eigenvalues
\label{sec:schur}


```{index} eigenvalue!conditioning of
```


```{index} condition number!of eigenvalues
```

Just as linear systems have condition numbers that quantify the effect of fixed precision, eigenvalue problems may be poorly conditioned too. While many possible results can be derived, we will use just one, the 
```{index} Bauer--Fike theorem
```
 **Bauer--Fike theorem**.\footnote{We will apply it only in the 2-norm, though it is more generally true.}  
::::proof:theorem}
  \label{thm:bauerfike}
  Let $\mathbf{A}\in\mathbb{C}^{n\times n}$ be diagonalizable, $\mathbf{A}=\mathbf{V}\mD\mathbf{V}^{-1}$, with
  eigenvalues $\lambda_1,\ldots,\lambda_n$. If $\mu$ is an eigenvalue
  of $\mathbf{A}+\mathbf{E}$ for a complex matrix $\mathbf{E}$, then
  
```{math}
    :label: bauerfike
    \min_{j=1,\ldots,n} |\mu - \lambda_j| \le \kappa(\mathbf{V}) \| \mathbf{E} \|,
  ```

  where $\|\cdot\|$ and $\kappa$ are in the 2-norm.
::::


% 
::::{proof:proof}
%   If the matrix $D-\mu I$ is singular, then $\mu=\lambda_j$ for some
%   $j$ and we are done. Otherwise, the eigenvalue condition on $\mu$
%   implies that $\mathbf{A}+\Delta A - \mu I$ is singular. We can multiply on
%   the left by $V^{-1}$ and on the right by $\mathbf{V}$ without changing the
%   singularity of the matrix. Using $\mathbf{A}=VDV^{-1}$, we have that
%   
```{math}
%     D + V^{-1}(\Delta A) V - \mu I
%   ```

%   is singular. So there exists a unit vector $\mathbf{x}$ such that
%   \begin{align*}
%     (D-\mu I)\mathbf{x} &= -V^{-1}(\Delta A)V \mathbf{x} \\
%     -(D-\mu I)^{-1} V^{-1}(\Delta A)V \mathbf{x} &= \mathbf{x}.
%   \end{align*}
%   By the definition of an induced matrix norm (see {eq}`normineq1`),
%   
```{math}
%     \| -(D-\mu I)^{-1 \| V^{-1}(\Delta A)V} \ge 1.
%   ```

%   Using norm inequality {eq}`normineq2`, we deduce
%   
```{math}
%     1 \le \| -(D-\mu I)^{-1 \|}\, \| V^{-1 \|} \,\| \Delta A \|\,
%     \| V \| = \| (D-\mu I)^{-1 \|}\, \| \Delta A \| \kappa(\mathbf{V}).
%   ```

%   The proof is finished by noting that, by diagonality (see
%   {ref}`prob-norms-diagnorm`),
%   
```{math}
%     \| (D-\mu I)^{-1 \|} = \max_{j} \frac{1}{|\lambda_j-\mu|}
%     = \frac{1}{\displaystyle \min_j |\lambda_j-\mu|}.
%   ```

% ::::


The Bauer--Fike theorem tells us that eigenvalues can be perturbed by an amount that is $\kappa(\mathbf{V})$ times larger than perturbations to the matrix. This result is a bit less straightforward than it might seem---eigenvectors are not unique, so there are multiple possible values for $\kappa(\mathbf{V})$. Still, the theorem indicates caution when a matrix has eigenvectors that form an ill-conditioned matrix. The limiting case of $\kappa(\mathbf{V})=\infty$ might be interpreted as indicating a nondiagonalizable matrix $\mathbf{A}$.

At the other extreme, if a unitary eigenvector matrix $\mathbf{V}$ can be found, then $\kappa(\mathbf{V})=1$ and {eq}`bauerfike` guarantees that eigenvalues are robust under perturbations to the original matrix $\mathbf{A}$. Such matrices are called 
```{index} matrix!normal
```
 **normal**, and they include the hermitian (or real symmetric) matrices. We consider them again in \secref{symm-eig}.


::::{proof:example}
  \inputexample{matrixanaly}{bauerfike}
::::


## Computing the EVD

In elementary linear algebra you use the characteristic polynomial to compute the eigenvalues of small matrices. However, computing polynomial roots in finite time is impossible for degree~5 and over. In principle one could use Newton-like methods to find all of the roots, but doing so is relatively slow and difficult. Furthermore we know that polynomial roots tend to become poorly conditioned when roots get close to one another (see {ref}`example-quadrootcond`).

Practical algorithms for computing the EVD go beyond the scope of this book. The essence of the matter is the connection to matrix powers, that is, $\mathbf{A}^k = \mathbf{V} \mathbf{D}^k \mathbf{V}^{-1}$. (We will see much more about the importance of matrix powers in \charef{krylov}.) If the eigenvalues have different complex magnitudes, then as $k\to\infty$ the entries on the diagonal of $\mathbf{D}^k$ become increasingly well separated and easy to pick out. It turns out that there is an astonishingly easy and elegant way to accomplish this separation without explicitly computing the matrix powers.


::::{proof:example}
  \inputexample{matrixanaly}{qriter}
::::


The process demonstrated in {ref}`example-qriter` is known as the 
```{index} Francis QR iteration
```
 *Francis QR iteration*, and it can be formulated as an $O(n^3)$ algorithm for finding the EVD. Such an algorithm is what the "eig" command uses.



% ## The Schur decomposition

% Mathematically, a matrix $\mathbf{V}$ whose columns are independent eigenvectors can be used to transform $\mathbf{A}$ to a diagonal matrix having the same eigenvectors: $V^{-1}AV=D$. The possibility of a large $\kappa(\mathbf{V})$, or even a nondiagonalizable matrixthe ill conditioning---or even nonexistence---of $\mathbf{V}$, however, suggests looking for alternative similarity transformations.

% A unitary matrix has a condition number equal to one, the lowest possible value. Therefore let us consider the idea of finding a unitary $\mathbf{U}$ such that $\mathbf{U}^{-1}\mathbf{A}\mU=\mathbf{U}^*\mathbf{A}\mU$ has easily computed eigenvalues. An important theorem tells us that this is always possible.
% 
::::proof:theorem}
%   \label{thm:schur}
%   For any complex square matrix $\mathbf{A}$ there exist a unitary $\mathbf{U}$ and an
%   upper triangular $T$ such that
%   
```{math}
%     :label: schur
%     \mathbf{A} = \mathbf{U}\m{T}\mathbf{U}^{-1} = \mathbf{U}\m{T}\mathbf{U}^*.
%   ```

%   That is, every square matrix is unitarily similar to a triangular matrix.
% ::::


% Equation {eq}`schur` is called a **Schur decomposition**. The value of the Schur decomposition stems from the fact that the eigenvalues of a triangular matrix $\m{T}$ are the diagonal entries of $\m{T}$. Hence unitary similarity transformations of $\mathbf{A}$ can be used to reveal its eigenvalues.

\begin{exercises}
  \input{matrixanaly/exercises/EigenvalueDecomposition}
\end{exercises}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\clearpage

# Singular value decomposition
\label{sec:svd}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

We now introduce another factorization that is as fundamental as the EVD.


::::proof:theorem}
  Let $\mathbf{A}\in\mathbb{C}^{m\times n}$. Then $\mathbf{A}$ can be written as
  
```{math}
    :label: svd
    \mathbf{A} = \mathbf{U} \mathbf{S} \mathbf{V}^*,
  ```

  where $\mathbf{U}\in\mathbb{C}^{m\times m}$ and $\mathbf{V}\in\mathbb{C}^{n\times n}$ are
  
```{index} unitary matrix
```
 unitary and $\mathbf{S}\in\mathbb{R}^{m\times n}$ is real and diagonal with
  nonnegative entries. If $\mathbf{A}$ is real, then so are $\mathbf{U}$ and
  $\mathbf{V}$ (which are then orthogonal matrices).
::::


\texthighlight{svd}{Equation {eq}`svd` is called a 
```{index} matrix!factorization|seealso {singular value decomposition
```
} 
```{index} singular value decomposition
```
 {term}`singular value decomposition`, or SVD, of $\mathbf{A}$.} The columns of $\mathbf{U}$ and $\mathbf{V}$ are called **left** and **right singular vectors**, respectively. The diagonal entries of $\mathbf{S}$, written $\sigma_1,\ldots,\sigma_r$, for $r=\min\{m,n\}$, are called the 
```{index} singular value
```
 **singular values** of $\mathbf{A}$. By convention the singular values are ordered so that

```{math}
  :label: svdorder
  \sigma_1 \ge \sigma_2 \ge \cdots \ge \sigma_r\ge 0, \qquad r=\min\{m,n\}.
```

We call $\sigma_1$ the 
```{index} singular value!principal
```
 **principal singular value** and $\mathbf{u}_{1}$ and $\mathbf{v}_{1}$ the **principal singular vectors**. The matrix $\mathbf{S}$ in the SVD is uniquely defined when the ordering is imposed, but the singular vectors are not---one could replace both $\mathbf{U}$ and $\mathbf{V}$ by their negatives, for example.


::::{proof:example}
  Suppose $\mathbf{A}$ is a real matrix and that $\mathbf{A}=\mathbf{U}\mS\mathbf{V}^T$ is an SVD. Then $\mathbf{A}^T=\mathbf{V}\mS^T\mathbf{U}^T$ meets all the requirements of an SVD for $\mathbf{A}^T$: the first and last matrices are orthogonal, and the middle matrix is diagonal with nonnegative entries. Hence $\mathbf{A}$ and $\mathbf{A}^T$ have the same singular values.
::::


## Interpreting the SVD

Another way to write $\mathbf{A}=\mathbf{U}\mS\mathbf{V}^*$ is as $\mathbf{A}\mV=\mathbf{U}\mS$. Taken columnwise, this means

```{math}
  :label: svdcolumns
  \mathbf{A} \mathbf{v}_{k} = \sigma_k \mathbf{u}_{k}, \qquad k=1,\ldots,r=\min\{m,n\}.
```

In words, each right singular vector is mapped by $\mathbf{A}$ to a scaled version of its corresponding left singular vector; the magnitude of scaling is its singular value. Together with the orthonormality of the columns of $\mathbf{U}$ and $\mathbf{V}$, these equations provide a simple graphical interpretation of the SVD. A representative $3\times 3$ case is shown in \figref{svdvisual}.

\begin{figure}[tbp]
  \centering
  \includegraphics[width=2in]{svdright}
  \hspace{8mm}
  \includegraphics[width=2in]{svdleft}
  \caption{Visual interpretation of the SVD for $3\times 3$ real matrices. Under multiplication by the matrix, a unit sphere, expressed in the right singular vector coordinates (shown at left), maps to an ellipsoid, given in the left singular vector coordinates (shown at right).}
  \label{fig:svdvisual}
\end{figure}

Both the SVD and the EVD describe a matrix in terms of some special vectors and a small number of scalars. {ref}`tab-evdsvd` summarizes the key differences. The SVD sacrifices having the same basis in both source and image spaces---after all, they may not even have the same dimension---but as a result gains orthogonality in both spaces.
\begin{table}[tb]
  \centering
  \caption{Differences between the EVD and the SVD.}
  \label{tab:evdsvd}
  \begin{tabular}{cc}
    **EVD**  & **SVD** \\\hline
    most square matrices & all rectangular and square matrices \\
    $\mathbf{A}\mathbf{x}_k = \lambda_k \mathbf{x}_k$    &   $\mathbf{A} \mathbf{v}_k = \sigma_k \mathbf{u}_k$ \\
    same basis for domain and range of $\mathbf{A}$ & two orthogonal bases \\
    may have poor conditioning & perfectly conditioned
  \end{tabular}
\end{table}




## SVD and the 2-norm

\texthighlight{svd2norm}{The SVD is intimately connected to the 2-norm,} as the following theorem describes. 
```{index} norm!matrix
```


::::proof:theorem}
  \label{thm:svd}
  Let $\mathbf{A}\in\mathbb{C}^{m\times n}$ have an SVD $\mathbf{A}=\mathbf{U}\mS\mathbf{V}^*$ in
  which {eq}`svdorder` holds. Then
  \begin{remunerate}
  \item The 2-norm satisfies
    
```{math}
      :label: svdnorm
      \| \mathbf{A} \|_2 = \sigma_1.
    ```

  \item The rank of $\mathbf{A}$ is the number of nonzero singular values.
  \item Let $r=\min\{m,n\}$. Then
    
```{math}
      :label: svdcond
      \kappa_2(\mathbf{A}) = \|\mathbf{A}\|_2\|\mathbf{A}^+\|_2 = \frac{\sigma_1}{\sigma_r},
    ```

    where a division by zero implies that $\mathbf{A}$ does not have full rank.
  % \item Let
  % 
```{math}
  %     :label: svdrank1
  %     \mA_k = \sigma_1 \mathbf{u}_{1}\mathbf{v}_{1}^* + \sigma_2
  %     \mathbf{u}_{2}\mathbf{v}_{2}^*
  %     + \cdots + \sigma_k \mathbf{u}_{k}\mathbf{v}_{k}^* = \mU_k\mS_k\mV_k^*,
  % ```

  % where $\mU_k$ and $\mV_k$ are $n\times k$ matrices with the first $k$
  % columns of $\mathbf{U}$ and $\mathbf{V}$, and $\mS_k$ is the upper-left $k\times k$
  % submatrix of $\mathbf{S}$. Then $\mA_n=\mathbf{A}$ and
  % 
```{math}
  %   \| \mathbf{A} - \mA_k \|_2 = \sigma_{k+1}, \qquad k=1,\ldots,n-1.
  % ```

  \end{remunerate}
::::


The conclusion {eq}`svdnorm` is intuitively clear from \figref{svdvisual}, and some straightforward vector calculus provides a proof (see {ref}`prob-svd-svdnormproof`). In the square case $m=n$, $\mathbf{A}$ having full rank is identical to being nonsingular.  The SVD is the usual means for computing the 2-norm and condition number of a matrix. We demonstrate this using MATLAB's {term}`svd` command.


::::{proof:example}
  \inputexample{matrixanaly}{svdproperties}
::::


## Connections to the EVD


```{index} matrix!factorization!eigenvalue
```

Let $\mathbf{A}=\mathbf{U}\mS\mathbf{V}^*$ be $m\times n$, and consider the square 
```{index} matrix!hermitian
```
 hermitian matrix $\mathbf{B}=\mathbf{A}^*\mathbf{A}$:

```{math}
  \mathbf{B} = (\mathbf{V}\mS^*\mathbf{U}^*) (\mathbf{U}\mS\mathbf{V}^*) = \mathbf{V}\mS^*\mathbf{S}\mV^* = \mathbf{V}(\mathbf{S}^T\mathbf{S})\mathbf{V}^{-1}.
```

Note that $\mathbf{S}^T\mathbf{S}$ is a diagonal $n \times n$ matrix. There are two cases to consider:
\begin{center}
  \begin{tabular}{lcl}
    \underline{$m \ge n$:} & $\qquad$ & \underline{$m < n$:} \\
    $\displaystyle \mathbf{S}^T\mathbf{S} =
    \begin{bmatrix}
      \sigma_1^2 & & \\
      & \ddots & \\
      & & \sigma_n^2
    \end{bmatrix}$, &  &
    $\displaystyle \mathbf{S}^T\mathbf{S} =
    \begin{bmatrix}
      \sigma_1^2 & & & \\
      & \ddots & & \\
      & & \sigma_m^2 & \\
      & & & \bm{0}
    \end{bmatrix}$,
  \end{tabular}
\end{center}
where the lower-right zero in the last matrix is $n-m$ square. In both cases we may conclude that the squares of the singular values of $\mathbf{A}$ are all eigenvalues of $\mathbf{B}$. Conversely, an EVD of $\mathbf{B}$ reveals the singular values and a set of right singular vectors of $\mathbf{A}$. The left singular vectors could then be deduced from the identity $\mathbf{A}\mV = \mathbf{U}\mS$.

Another close connection between EVD and SVD comes via the $(m+n)\times (m+n)$ matrix

```{math}
  :label: svdaugment
  \mathbf{C} =
  \begin{bmatrix}
    0 & \mathbf{A}^* \\ \mathbf{A} & 0
  \end{bmatrix}.
```

If $\sigma$ is a singular value of $\mathbf{B}$, then $\sigma$ and $-\sigma$ are eigenvalues of $\mathbf{C}$, and the associated eigenvector immediately reveals a left and a right singular vector (see {ref}`prob-svd-svdtoevd`). This connection is implicitly exploited by software to compute the SVD.

## Thin form

In \secref{QR} we saw that a matrix has both a "full" and a "thin" or reduced form of the QR factorization. 
```{index} singular value decomposition!thin form
```
 A similar situation holds with the SVD. Suppose $\mathbf{A}$ is $m\times n$ with $m > n$ and let $\mathbf{A}=\mathbf{U}\mS\mathbf{V}^*$ be an SVD. The last $m-n$ rows of $\mathbf{S}$ are all zero due to the fact that $\mathbf{S}$ is diagonal. Hence
\begin{align}
  \mathbf{U} \mathbf{S} & =
  \begin{bmatrix}
    \mathbf{u}_1 & \cdots & \mathbf{u}_n & \mathbf{u}_{n+1} & \cdots & \mathbf{u}_m
  \end{bmatrix}
  \begin{bmatrix}
    \sigma_1 & &  \\
    & \ddots &  \\
    & & \sigma_n \\
    & & \\
    & \bm{0} & \\
    & &
  \end{bmatrix} \notag\\
  &=
  \begin{bmatrix}
    \mathbf{u}_1 & \cdots & \mathbf{u}_n
  \end{bmatrix}
  \begin{bmatrix}
    \sigma_1 & &  \\
    & \ddots &  \\
    & & \sigma_n
  \end{bmatrix} = \widehat{\mathbf{U}} \widehat{\mathbf{S}},
        :label: svd-reduced
\end{align}
in which $\widehat{\mathbf{U}}$ is $m\times n$ and $\widehat{\mathbf{S}}$ is $n\times n$. This allows us to define the **thin SVD** $\mathbf{A}=\widehat{\mathbf{U}}\widehat{\mathbf{S}}\mathbf{V}^*$, in which $\widehat{\mathbf{S}}$ is square and diagonal and $\widehat{\mathbf{U}}$ is 
```{index} ONC
```
 ONC but not unitary. This form is computationally preferable when $m\gg n$, since it requires far less storage and contains the same information about $\mathbf{A}$. In MATLAB one uses \verb|svd(A,0)| to get the thin form. In the dual case where $m<n$, you can use the reduced SVD of $\mathbf{A}^*$ and then take the hermitian to get a reduced SVD for $\mathbf{A}$.


\begin{exercises}
	\input{matrixanaly/exercises/SVD}
\end{exercises}



%\clearpage
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% # Computing the decompositions
% \label{sec:compdecomp}
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %(This section can be skipped without interrupting the flow or breaking
% %dependencies.)


% There are some other important ways to connect singular values and eigenvalues.

\clearpage

# Symmetry and definiteness
\label{sec:symm-eig}


```{index} matrix!symmetric
```

As we saw in \secref{specialmatrix}, symmetry can simplify the LU factorization into a symmetric form, $\mathbf{A}=\mathbf{L}\mD\mathbf{L}^T$. Certain specializations occur too for the eigenvalue and singular value factorizations. In this section we stay with complex-valued matrices, so we are interested in the case when $\mathbf{A}^*=\mathbf{A}$, or $\mathbf{A}$ is hermitian. However, we often loosely speak of "symmetry" to mean this property in the complex case. All of the statements in this section easily specialize to the real case.

## Unitary diagonalization

Suppose now that $\mathbf{A}^*=\mathbf{A}$ and that $\mathbf{A}=\mathbf{U}\mS\mathbf{V}^*$ is an SVD. Since $\mathbf{S}$ is real and square, we have

```{math}
  \mathbf{A}^* = \mathbf{V} \mathbf{S}^* \mathbf{U}^* = \mathbf{V} \mathbf{S} \mathbf{U}^*,
```

and it's tempting to conclude that $\mathbf{U}=\mathbf{V}$. Happily, this is nearly true. The following theorem is typically proved in an advanced linear algebra course.

::::proof:theorem}[Spectral decomposition]
  \label{thm:spec-decomp} If $\mathbf{A}=\mathbf{A}^*$, then $\mathbf{A}$ has a diagonalization $\mathbf{A}=\mathbf{V} \mathbf{D} \mathbf{V}^{-1}$ in which $\mathbf{V}$ is unitary and $\mathbf{D}$ is diagonal and real.
::::

Another way to state the result of \thmref{spec-decomp} is that \texthighlight{hermitianevec}{a hermitian matrix has a complete set of orthonormal eigenvectors---that is, a *unitary diagonalization*---and real eigenvalues.} In this case, the EVD $\mathbf{A}=\mathbf{V}\mD\mathbf{V}^{-1}=\mathbf{V} \mathbf{D} \mathbf{V}^*$ is almost an SVD.

::::proof:theorem}
  If $\mathbf{A}^*=\mathbf{A}$ and $\mathbf{A}=\mathbf{V}\mD\mathbf{V}^{-1}$ is a 
```{index} unitary matrix
```
 unitary diagonalization, then
  
```{math}
    :label: herm-svd
    \mathbf{A} = (\mathbf{V}\m{T})\cdot |\mathbf{D}|\cdot \mathbf{V}^*
  ```

  is an SVD, where $|\mathbf{D}|$ is the elementwise absolute value and $\m{T}$ is diagonal with $|T_{ii}|=1$ for all $i$.
::::


::::{proof:proof}
  Let $T_{ii}=\operatorname{sign}(D_{ii})$. Then $\m{T}^2=\mathbf{I}$ and $|\mathbf{D}|=\m{T}\mathbf{D}$. The result follows.
::::


The converse of \thmref{spec-decomp} is also true: every matrix with a unitary diagonalization and real eigenvalues is hermitian. However, there are nonhermitian matrices that meet just the requirement of a unitary EVD; any such matrix is called **normal**.


::::{proof:example}
  \inputexample{matrixanaly}{svdnormal}
::::


Now consider again \thmref{bauerfike}, which says that the condition number of the eigenvalues is bounded above by $\kappa(\mathbf{V})$, where $\mathbf{V}$ is an eigenvector matrix. Because $\kappa=1$ for any unitary or orthogonal matrix, \thmref{spec-decomp} then implies that the condition number of the eigenvalues of a hermitian or any normal matrix is one. That is, \texthighlight{normaleval}{eigenvalues of a normal matrix can be changed by no more than the norm of the perturbation to the matrix.}


::::{proof:example}
  \inputexample{matrixanaly}{normalperturb}
::::


## Rayleigh quotient


```{index} Rayleigh quotient
```

Recall that for a matrix $\mathbf{A}$ and compatible vector $\mathbf{x}$, the quadratic form $\mathbf{x}^* \mathbf{A} \mathbf{x}$ is a scalar. With a suitable normalization, it becomes the {term}`Rayleigh quotient`

```{math}
  :label: rayleigh
  R_{\mathbf{A}}(\mathbf{x}) = \frac{ \mathbf{x}^* \mathbf{A} \mathbf{x}}{\mathbf{x}^* \mathbf{x}}.
```

If $\mathbf{v}$ is an eigenvector such that $\mathbf{A} \mathbf{v}=\lambda \mathbf{v}$, then one easily calculates that $R_{\mathbf{A}}(\mathbf{v})=\lambda.$ That is, \texthighlight{rayleighquotient}{the Rayleigh quotient maps an eigenvector into its associated eigenvalue.}

If $\mathbf{A}^*=\mathbf{A}$, then the Rayleigh quotient has another interesting property: $\nabla R_{\mathbf{A}}(\mathbf{v})=oldsymbol{0}$ if $\mathbf{v}$ is an eigenvector. By a multidimensional Taylor series, then,

```{math}
  :label: rq-series
  R_{\mathbf{A}}(\mathbf{v}+\epsilon\mathbf{z}) = R_{\mathbf{A}}(\mathbf{v}) + 0 + O( \epsilon^2) =  \lambda + O( \epsilon^2),
```

as $\epsilon\to 0$. The conclusion is that a good estimate of an eigenvector becomes an even better estimate of an eigenvalue.


::::{proof:example}
  \inputexample{matrixanaly}{rayquo}
::::


## Definite and indefinite matrices


```{index} matrix!positive definite
```

In the real case, we called a symmetric matrix $\mathbf{A}$ *symmetric positive definite* (SPD)
if $\mathbf{x}^T \mathbf{A}\mathbf{x} > 0 $ for all nonzero vectors $\mathbf{x}$. In the complex case the relevant property is 
```{index} hermitian positive definite
```
 **hermitian positive definite** (HPD), meaning that $\mathbf{A}^*=\mathbf{A}$ and $\mathbf{x}^* \mathbf{A}\mathbf{x} > 0$ for all complex vectors $\mathbf{x}$. Putting this property together with the Rayleigh quotient leads to

::::proof:theorem}
  \label{thm:hpd} If $\mathbf{A}^*=\mathbf{A}$, then the following statements are equivalent.
  \begin{remunerate}
  \item $\mathbf{A}$ is HPD.
  \item The eigenvalues of $\mathbf{A}$ are positive numbers.
  \item Any unitary EVD of $\mathbf{A}$ is also an SVD of $\mathbf{A}$.
  \end{remunerate}
::::


Naturally, a hermitian matrix with all negative eigenvalues is called *negative definite*, and one with eigenvalues of different signs is *indefinite*. Finally, if one or more eigenvalues is zero and the rest have one sign, it is *semidefinite*.


\begin{exercises}
  \input{matrixanaly/exercises/Symmetry}
\end{exercises}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\clearpage
# Dimension reduction
\label{sec:dimreduce}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


```{index} singular value decomposition
```


```{index} dimension reduction
```

The SVD has another important property that proves very useful in a variety of applications. Let $\mathbf{A}$ be a real $m\times n$ matrix with SVD $\mathbf{A}=\mathbf{U}\mS\mathbf{V}^T$ and (momentarily) $m\ge n$. Another way of writing the thin form of the SVD is
\begin{align}
  \mathbf{A} = \widehat{\mathbf{U}}\widehat{\mathbf{S}}\mathbf{V}^T &=
  \begin{bmatrix}
    \mathbf{u}_1 & \mathbf{u}_2 & \cdots & \mathbf{u}_n
  \end{bmatrix}
  \begin{bmatrix}
    \sigma_1 & & \\
    & \ddots & \\
    & & \sigma_n
  \end{bmatrix}
        \begin{bmatrix}
          \mathbf{v}_1^T \\ \vdots \\ \mathbf{v}_n^T
        \end{bmatrix}  \notag \\
  &=
  \begin{bmatrix}
    \sigma_1\mathbf{u}_1  & \cdots & \sigma_n\mathbf{u}_n
  \end{bmatrix}
                               \begin{bmatrix}
                                  \mathbf{v}_1^T \\ \vdots \\ \mathbf{v}_n^T
                               \end{bmatrix} \notag \\
  :label: svdsum
  &= \sigma_1 \mathbf{u}_{1}\mathbf{v}_{1}^T + \cdots + \sigma_r \mathbf{u}_{r}\mathbf{v}_{r}^T = \sum_{i=1}^r \sigma_i \mathbf{u}_{i}\mathbf{v}_{i}^T,
\end{align}
where $r$ is the rank of $\mathbf{A}$. The final formula also holds for the case $m<n$.

Each outer product $\mathbf{u}_{i}\mathbf{v}_{i}^T$ is a rank-1 matrix of unit norm. Thanks to the ordering of singular values, then, equation {eq}`svdsum` expresses $\mathbf{A}$ as a sum of decreasingly important contributions. This motivates the definition, for $1\le k \le r$,

```{math}
  :label: svdlowrank
  \mA_k = \sum_{i=1}^k \sigma_i \mathbf{u}_{i}\mathbf{v}_{i}^T = \mU_k \mS_k \mV_k^T.
```

where $\mU_k$ and $\mV_k$ are the first $k$
columns of $\mathbf{U}$ and $\mathbf{V}$, respectively, and $\mS_k$ is the upper-left $k\times k$
submatrix of $\mathbf{S}$.

The rank of a sum of matrices is always less than or equal to the sum of the ranks, so $\mA_k$ is a rank-$k$ approximation to $\mathbf{A}$. It turns out that \texthighlight{bestrank}{$\mA_k$ is the *best* rank-$k$ approximation of $\mathbf{A}$,} as measured in the matrix 2-norm.

::::proof:theorem}
  \label{thm:best-rank-k}
  Suppose $\mathbf{A}$ has rank $r$ and let $\mathbf{A}=\mathbf{U}\mS\mathbf{V}^T$ be an SVD. Let $\mA_k$ be as in {eq}`svdlowrank` for $1\le k < r$. Then
  \begin{enumerate}
  \item $\| \mathbf{A} - \mA_k \|_2 = \sigma_{k+1}, \qquad k=1,\ldots,r-1$.
  \item If the rank of $\mathbf{B}$ is $k$ or less, then $\| \mathbf{A}-\mathbf{B} \|_2\ge \sigma_{k+1}$.
  \end{enumerate}
::::


::::{proof:proof}[Proof of part 1]
  Note that {eq}`svdlowrank` is identical to {eq}`svdsum` with $\sigma_{k+1},\ldots,\sigma_r$ all set to zero. This implies that
  
```{math}
    \mathbf{A} - \mA_k = \mathbf{U}(\mathbf{S}-\hat{\mathbf{S}})\mathbf{V}^T,
  ```

  where $\hat{\mathbf{S}}$ has those same values of $\sigma_i$ replaced by zero. But that makes the above an SVD of $\mathbf{A} - \mA_k$, with singular values $0,\ldots,0,\sigma_{k+1},\ldots,\sigma_r$, the largest of which is $\sigma_{k+1}$. That proves the first claim.
::::



%
%\autoref{fig:dimreduce} shows a geometric interpretation: the vectors $\sigma_i \mathbf{u}_i$ are semi-axes of an $r$-dimensional hyperellipse that forms the image of the unit hypersphere under the mapping $\mathbf{x}\mapsto\mathbf{A}\mathbf{x}$, and the low-rank approximation~\autoref{eq:svdlowrank} sets the shortest $r-k$ semi-axes to zero, effectively flattening out those dimensions.
% \begin{figure}
%   \centering
%   \includegraphics[width=3.5in]{dimreduce}
%   \caption{Dimension reduction via the SVD. The left singular vectors give the principal directions of the ellipsoid that results from mapping the unit sphere. Zeroing out the direction with the smallest singular value leads to the ellipse in the figure, lying in the plane shown. Zeroing out the next smallest direction leads to a line in space, which is the best possible one-dimensional approximation of the range of this matrix.}
%   \label{fig:dimreduce}
% \end{figure}
%
If the singular values of $\mathbf{A}$ decrease sufficiently rapidly, then $\mA_{k}$ may capture the "most significant" behavior of the matrix for a reasonably small value of $k$.


```{index} matrix!as image
```


::::{proof:example}
  \inputexample{matrixanaly}{hellosvd}
::::


## Capturing major trends

The use of dimension reduction offered by low-rank SVD approximation goes well beyond simply reducing computation time. By isolating the most important contributions to the matrix, \texthighlight{dimreduce}{dimension reduction can uncover deep connections and trends that are otherwise obscured by lower-order effects and noise.}

One useful way to quantify the decay in the singular values is to compute

```{math}
  :label: sing-val-decay
  \tau_k = \frac{\sum_{i=1}^k \sigma_i^2}{\sum_{i=1}^r \sigma_i^2}, \quad k=1,\ldots,r.
```

Clearly $0\le \tau_k \le 1$ and $\tau_k$ is nondecreasing as a function of $k$. We can think of $\tau_k$ as the fraction of "energy" contained in the singular values up to and including the $k$th.\footnote{In statistics this quantity may be interpreted as the fraction of explained variance.}


::::{proof:example}
  \inputexample{matrixanaly}{voting}
::::


Not all data sets can be reduced effectively to a small number of dimensions, but as {ref}`example-voting` shows, in some cases reduction reveals information that may correspond to real-world understanding.

\begin{exercises}
	\input{matrixanaly/exercises/DimReduce}
\end{exercises}

% 
::::{proof:example}
%   \label{exa:mandrills}
%   Here we show how SVD compression works on an image that ships with
%   MATLAB. In this case the image is stored as a matrix of integer
%   values that we interpret as gray levels.
% \begin{verbatim}
% >> load mandrill
% >> size(X)

% ans =
%    480   500
% >> [U,S,V] = svd(X);
% >> X1 = U(:,1:40)*S(1:40,1:40)*V(:,1:40)';
% >> colormap(gray(220))
% >> subplot(121), image(X)
% >> subplot(122), image(X1)
% \end{verbatim}
%   \begin{center}
%     \includegraphics{mandrills}
%   \end{center}
%   SVD compression isn't competitive with standards like JPEG for
%   single images in terms of preserving fine detail, but it functions
%   quite well in some systems for facial recognition.
% ::::


\clearpage
\subsection*{Key ideas in this chapter}
\begin{remunerate}
\item Tables, graphs, and images are all examples of two-dimensional data that can be usefully interpreted as matrices (\hiref{matrixinsight}).
\item A square matrix may have an eigenvalue decomposition $\mathbf{A}=\mathbf{V}\mD\mathbf{V}^{-1}$, where $\mathbf{D}$ is diagonal (\hiref{evd}).
\item Multiplication by the inverse of a matrix performs a change of basis into the coordinates associated with the matrix. An EVD changes basis to one in which the matrix acts diagonally (\hiref{changebasis}).
\item Any matrix has a singular value decomposition $\mathbf{A}=\mathbf{U}\mS\mathbf{V}^*$, where $\mathbf{S}$ is diagonal and nonnegative, and $\mathbf{U}$ and $\mathbf{V}$ are unitary (\hiref{svd}).
\item The SVD is intimately connected to the 2-norm (\hiref{svd2norm}).
\item A hermitian matrix has a complete set of orthonormal eigenvectors and real eigenvalues (\hiref{hermitianevec}).
\item The eigenvalues of a normal matrix are perturbed by no more than the size of the perturbation to the matrix (\hiref{normaleval}).
\item The Rayleigh quotient is a function mapping an eigenvector into its associated eigenvalue, and an eigenvector estimate into an eigenvalue estimate (\hiref{rayleighquotient}).
\item Truncating the SVD provides the best low-rank approximation to a matrix (\hiref{bestrank}).
\item Dimension reduction via the SVD can uncover connections that are otherwise obscured by lower-order effects and noise (\hiref{dimreduce}).
\end{remunerate}


\subsection*{Where to learn more}

Details on the computation of the eigenvalue and singular value decompositions are presented at length in~\cite{StewartVol2} and more briefly in Chapters~7 and~8 of~\cite{GolubVan96}. A classic reference on the particulars of the symmetric case is~\cite{Parlett1980}, while~\cite{TrefEmb05} focuses on the non-normal case. Dimension reduction via the SVD often goes by the name *principal component analysis*, which is the subject of~\cite{Jolliffe2002}.



%%% Local Variables:
%%% mode: latex
%%% TeX-master: "natext"
%%% End:
