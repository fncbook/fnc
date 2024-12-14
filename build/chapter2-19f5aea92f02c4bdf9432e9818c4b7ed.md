---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.4
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

```{code-cell} ipython3
# from scipy import *
import numpy as np
from matplotlib.pyplot import *
from scipy.linalg import solve
# from numpy.linalg import *
import scipy.sparse as sparse
from scipy.sparse.linalg import splu
from timeit import default_timer as timer
import FNC
```

```{code-cell} ipython3
# This (optional) block is for improving the display of plots.
# from IPython.display import set_matplotlib_formats
# set_matplotlib_formats("svg","pdf")
# %config InlineBackend.figure_format = 'svg'
rcParams["figure.figsize"] = [7, 4]
rcParams["lines.linewidth"] = 2
rcParams["lines.markersize"] = 4
rcParams['animation.html'] = "jshtml"  # or try "html5"
```

<!-- SECTION 1 -->
(demo-interp-vander-python)=
``````{dropdown} Linear system for polynomial interpolation
:open: false

We create two vectors for data about the population of China. The first has the years of census data and the other has the population, in millions of people.

```{code-cell} ipython3
year = np.arange(1980, 2020, 10)   # from 1980 to 2020 by 10
pop = np.array([984.736, 1148.364, 1263.638, 1330.141])
```

It's convenient to measure time in years since 1980. 

```{code-cell} ipython3
t = year - 1980
y = pop
```

Now we have four data points $(t_1,y_1),\dots,(t_4,y_4)$, so $n=4$ and we seek an interpolating cubic polynomial. We construct the associated Vandermonde matrix: 

```{code-cell} ipython3
V = np.vander(t)
print(V)
```

To solve a linear system $\mathbf{V} \mathbf{c} = \mathbf{y}$ for the vector of polynomial coefficients, we use `solve` (imported from `scipy.linalg`):

```{code-cell} ipython3
c = solve(V, y)
print(c)
```
::::{grid} 1 1 2 2
The algorithms used by `solve` are the main topic of this chapter. As a check on the solution, we can compute the *residual* $\mathbf{y} - \mathbf{V} \mathbf{c}$, which should be small (near machine precision).

:::{card}
Matrix multiplication in NumPy is done with `@` or `matmul`.
:::
::::

```{code-cell} ipython3
print(y - V @ c)
```

By our definitions, the coefficients in `c` are given in descending order of power in $t$. We can use the resulting polynomial to estimate the population of China in 2005:

```{code-cell} ipython3
p = np.poly1d(c)          # construct a polynomial
print(p(2005 - 1980))     # apply the 1980 time shift
```

The official figure was 1303.72, so our result is rather good.

We can visualize the interpolation process. First, we plot the data as points. Then we add a plot of the interpolant, taking care to shift the $t$ variable back to actual years.

```{code-cell} ipython3
scatter(year, y, color="k", label="data");
tt = np.linspace(0, 30, 300)   # 300 times from 1980 to 2010
plot(1980 + tt, p(tt), label="interpolant");
xlabel("year");
ylabel("population (millions)");
title("Population of China");
legend();
```
``````

(demo-matrices-python)=
``````{dropdown} Matrix operations
:open: false

```{note}
While NumPy does have distinct representations for matrices and 2D arrays, use of the explicit matrix class is officially discouraged. We follow this advice here and use arrays to represent both matrices and vectors.
```

:::{index} ! Python; array, ! Python; shape
:::

<!-- :::{index}
see: Python; size, Python; shape
::: -->

A vector is created using square brackets and commas to enclose and separate its entries.

```{code-cell} ipython3
x = np.array([3, 3, 0, 1, 0 ])
print(x.shape)
```

To construct a matrix, you nest the brackets to create a "vector of vectors". The inner vectors are the rows.

```{code-cell} ipython3
A = np.array([ 
    [1, 2, 3, 4, 5],
    [50, 40, 30, 20, 10], 
    [np.pi, np.sqrt(2), np.exp(1), (1+np.sqrt(5))/2, np.log(3)] 
    ])

print(A)
print(A.shape)
```

In this text, we treat all vectors as equivalent to matrices with a single column. That isn't true in NumPy, because even an $n \times 1$ array has two dimensions, unlike a vector.

```{code-cell} ipython3
np.array([[3], [1], [2]]).shape
```

:::{index} ! Python; hstack, ! Python; vstack
:::

You can concatenate arrays with compatible dimensions using `hstack` and `vstack`.

```{code-cell} ipython3
print( np.hstack([A, A]) )
```

```{code-cell} ipython3
print( np.vstack([A, A]) )
```

```{index} ! Python; transpose, ! Python; adjoint
```

Transposing a matrix is done by appending `.T` to it. 

```{code-cell} ipython3
print(A.T)
```

For matrices with complex values, we usually want instead the adjoint or hermitian, which is `.conj().T`. 

```{code-cell} ipython3
print((x + 1j).conj().T)
```

:::{index} ! Python; arange, ! Python; linspace
:::

There are many convenient shorthand ways of building vectors and matrices other than entering all of their entries directly or in a loop. To get a vector with evenly spaced entries between two endpoints, you have two options.

```{code-cell} ipython3
print(np.arange(1, 7, 2))   # from 1 to 7 (not inclusive), step by 2        
```

```{code-cell} ipython3
print(np.linspace(-1, 1, 5))   # from -1 to 1 (inclusive), with 5 total values
```

The practical difference between these is whether you want to specify the step size in `arange` or the number of points in `linspace`.

Accessing an element is done by giving one (for a vector) or two index values in square brackets. **In Python, indexing always starts with zero, not 1.**

```{code-cell} ipython3
A = np.array([ 
    [1, 2, 3, 4, 5],
    [50, 40, 30, 20, 10], 
    np.linspace(-5, 5, 5) 
    ])
x = np.array([3, 2, 0, 1, -1 ])
```

```{code-cell} ipython3
print("row 2, col 3 of A:", A[1, 2])
print("first element of x:", x[0])
```

```{index} ! Python; slice, ! Python; \:
```
:::{index} ! Python; indexing arrays
:::

The indices can be ranges, in which case a **slice** or block of the matrix is accessed. You build these using a colon in the form `start:stop`. However, the last value of this range is `stop-1`, not `stop`.

```{code-cell} ipython3
print(A[1:3, 0:2])    # rows 2 and 3, cols 1 and 2
```

If `start` or `stop` is omitted, the range extends to the first or last index.

```{code-cell} ipython3
print(x[1:])  # elements 2 through the end
```

```{code-cell} ipython3
print(A[:2, 0])  # first two rows in column 1
```

Notice in the last case above that even when the slice is in the shape of a column vector, the result is just a vector with one dimension and neither row nor column shape.

There are more variations on the colon ranges. A negative value means to count from the end rather than the beginning. And a colon by itself means to include everything from the relevant dimension.

```{code-cell} ipython3
print(A[:-1, :])    # all rows up to the last, all columns
```

Finally, `start:stop:step` means to step size or stride other than one. You can mix this with the other variations.

```{code-cell} ipython3
print(x[::2])  # all the odd indexes
```

```{code-cell} ipython3
print(A[:, ::-1])  # reverse the columns
```

The matrix and vector senses of addition, subtraction, and scalar multiplication and division are all handled by the usual symbols. Two matrices of the same size (what NumPy calls shape) are operated on elementwise. 

```{code-cell} ipython3
print(A - 2 * np.ones([3, 5]))  # subtract two from each element
```

```{index} ! Python; broadcast
```

If one operand has a smaller number of dimensions than the other, Python tries to **broadcast** it in the "missing" dimension(s), and the operation proceeds if the resulting shapes are identical. 

```{code-cell} ipython3
print(A - 2)    # subtract two from each element
```

```{code-cell} ipython3
u = np.array([1, 2, 3, 4, 5])
print(A - u)    # repeat this row for every row of A
```

```{code-cell} ipython3
:tags: raises-exception
v = np.array([1, 2, 3])
print(A - v)  # broadcasting this would be 3x3, so it's an error
```

```{code-cell} ipython3
print(A - v.reshape([3, 1]))    # broadcasts to each column of A
```

```{index} ! Python; \@, ! Python; matmul
```

<!-- ```{index} 
see: Python; matrix multiplication, Python; \@
``` -->

```{index} ! Python; diag
```

Matrix–matrix and matrix–vector products are computed using `@` or `matmul`.

```{code-cell} ipython3
B = np.diag([-1, 0, -5])    # create a diagonal 3x3
print(B @ A)    # matrix product
```

$AB$ is undefined for these matrix sizes. 

```{code-cell} ipython3
:tags: raises-exception
print(A @ B)    # incompatible sizes
```

```{index} ! Python; elementwise multiplication, Python; broadcasting
```

The multiplication operator `*` is reserved for elementwise multiplication. Both operands have to be the same size, after any potential broadcasts.

```{code-cell} ipython3
:tags: raises-exception
print(B * A)    # not the same size, so it's an error
```

```{code-cell} ipython3
print((A / 2) * A)    # elementwise
```
To raise to a power elementwise, use a double star. This will broadcast as well.

```{code-cell} ipython3
print(B)
print(B**3)
```

```{code-cell} ipython3
print(x)
print(2.0**x)
```

```{danger}
If `A` is a matrix, `A**2` is *not* the same as mathematically raising it to the power 2.
```


```{index} Python; broadcasting
```

Most of the mathematical functions, such as cos, sin, log, exp and sqrt, expecting scalars as operands will be broadcast to arrays.

```{code-cell} ipython3
print(np.cos(np.pi * x))      
```
``````

<!-- SECTION 3 -->
(demo-systems-backslash-python)=
``````{dropdown} Solving linear systems
:open: false

For a square matrix $A$, the command `solve(A, B)` is mathematically equivalent to $\mathbf{A}^{-1} \mathbf{b}$. 

```{code-cell} ipython3
A = np.array([[1, 0, -1], [2, 2, 1], [-1, -3, 0]])
b = np.array([1, 2, 3])
```

```{code-cell} ipython3
x = solve(A, b)
print(x)
```

```{index} residual
```

One way to check the answer is to compute a quantity known as the **residual**. It is (ideally) close to machine precision(relative to the elements in the data). 

```{code-cell} ipython3
residual = b - A @ x
print(residual)
```

If the matrix $\mathbf{A}$ is singular, you may get an error.

```{code-cell} ipython3
:tags: raises-exception
A = np.array([[0, 1], [0, 0]])
b = np.array([1, -1])
solve(A, b)    # error, singular matrix
```

A linear system with a singular matrix might have no solution or infinitely many solutions, but in either case, a numerical solution becomes trickier. Detecting singularity is a lot like checking whether two floating-point numbers are *exactly* equal: because of roundoff, it could be missed. We're headed toward a more robust way to fully describe this situation.
``````

(function-forwardsub-python)=
``````{dropdown} Forward substitution
:open:
```{literalinclude} ../python/pkg/FNC/FNC02.py
:filename: forwardsub.py
:start-line: 2
:end-line: 14
:linenos: true
:language: python
```
``````

(function-backsub-python)=
``````{dropdown} Backward substitution
:open:
```{literalinclude} ../python/pkg/FNC/FNC02.py
:filename: backsub.py
:start-line: 17
:end-line: 29
:linenos: true
:language: python
```
``````

(demo-systems-triangular-python)=
``````{dropdown} Triangular systems of equations
:open: false

```{index} ! Python; tril, ! Python; triu
```

It's easy to get just the lower triangular part of any matrix using the `tril` function.

```{code-cell} ipython3
A = 1 + np.floor(9 * np.random.rand(5, 5))
L = np.tril(A)
print(L)
```

We'll set up and solve a linear system with this matrix.

```{code-cell} ipython3
b = np.ones(5)
x = FNC.forwardsub(L, b)
print(x)
```

It's not clear how accurate this answer is. However, the residual should be zero or comparable to $\macheps$.

```{code-cell} ipython3
b - L @ x
```

Next we'll engineer a problem to which we know the exact answer. 

```{code-cell} ipython3
alpha = 0.3;
beta = 2.2;
U = np.diag(np.ones(5)) + np.diag([-1, -1, -1, -1], k=1)
U[0, 3:5] = [ alpha - beta, beta ]
print(U)
```

```{code-cell} ipython3
x_exact = np.ones(5)
b = np.array([alpha, 0, 0, 0, 1])
x = FNC.backsub(U, b)
print("error:", x - x_exact)
```

Everything seems OK here. But another example, with a different value for $\beta$, is more troubling.

```{code-cell} ipython3
alpha = 0.3;
beta = 1e12;
U = np.diag(np.ones(5)) + np.diag([-1, -1, -1, -1], k=1)
U[0, 3:5] = [ alpha - beta, beta ]
b = np.array([alpha, 0, 0, 0, 1])

x = FNC.backsub(U, b)
print("error:", x - x_exact)
```

It's not so good to get 4 digits of accuracy after starting with sixteen! But the source of the error is not hard to track down. Solving for $x_1$ performs $(\alpha-\beta)+\beta$ in the first row. Since $|\alpha|$ is so much smaller than $|\beta|$, this a recipe for losing digits to subtractive cancellation.
``````

<!-- SECTION 4 -->
(demo-lu-outertri-python)= 
``````{dropdown} Triangular outer products
:open: false

```{index} Python; tril, Python; triu
```
We explore the outer product formula for two random triangular matrices.

```{code-cell}
from numpy.random import randint
L = np.tril(randint(1, 10, size=(3, 3)))
print(L)
```

```{code-cell}
U = np.triu(randint(1, 10, size=(3, 3)))
print(U)
```

Here are the three outer products appearing in the sum in {eq}`matrixouter`:

```{code-cell}
print(np.outer(L[:, 0], U[0, :]))
```

```{code-cell}
print(np.outer(L[:, 1], U[1, :]))
```

```{code-cell}
print(np.outer(L[:, 2], U[2, :]))
```

Simply because of the triangular zero structures, only the first outer product contributes to the first row and first column of the entire product. 
``````


(demo-lu-derive-python)=
``````{dropdown} LU factorization
:open:

For illustration, we work on a $4 \times 4$ matrix. We name it with a subscript in preparation for what comes.

```{code-cell} ipython3
A_1 = np.array([
     [2,    0,    4,    3], 
     [-4,    5,   -7,  -10], 
     [1,   15,    2,   -4.5],
     [-2,    0,    2,  -13]
        ])
L = np.eye(4)
U = np.zeros((4, 4));
```

Now we appeal to {eq}`outer-row1`. Since $L_{11}=1$, we see that the first row of $\mathbf{U}$ is just the first row of $\mathbf{A}_1$.

```{code-cell}
U[0, :] = A_1[0, :]
print(U)
```

From {eq}`outer-col1`, we see that we can find the first column of $\mathbf{L}$ from the first column of $\mathbf{A}_1$. 

```{code-cell}
L[:, 0] = A_1[:, 0] / U[0, 0]
print(L)
```

We have obtained the first term in the sum {eq}`matrixouter` for $\mathbf{L}\mathbf{U}$, and we subtract it away from $\mathbf{A}_1$.

```{code-cell}
A_2 = A_1 - np.outer(L[:, 0],  U[0, :])
```

Now $\mathbf{A}_2 = \boldsymbol{\ell}_2\mathbf{u}_2^T + \boldsymbol{\ell}_3\mathbf{u}_3^T + \boldsymbol{\ell}_4\mathbf{u}_4^T.$ If we ignore the first row and first column of the matrices in this equation, then in what remains we are in the same situation as at the start. Specifically, only $\boldsymbol{\ell}_2\mathbf{u}_2^T$ has any effect on the second row and column, so we can deduce them now.

```{code-cell}
U[1, :] = A_2[1, :]
L[:, 1] = A_2[:, 1] / U[1, 1]
print(L)
```
If we subtract off the latest outer product, we have a matrix that is zero in the first *two* rows and columns. 

```{code-cell}
A_3 = A_2 - np.outer(L[:, 1], U[1, :])
```

Now we can deal with the lower right $2\times 2$ submatrix of the remainder in a similar fashion.

```{code-cell}
U[2, :] = A_3[2, :]
L[:, 2] = A_3[:, 2] / U[2, 2]
A_4 = A_3 - np.outer(L[:, 2], U[2, :])
```

Finally, we pick up the last unknown in the factors.

```{code-cell}
U[3, 3] = A_4[3, 3]
```

We now have all of $\mathbf{L}$,

```{code-cell} 
print(L)
```

and all of $\mathbf{U}$,

```{code-cell}
print(U)
```

We can verify that we have a correct factorization of the original matrix by computing the backward error:

```{code-cell} 
A_1 - L @ U
```

In floating point, we cannot expect the difference to be exactly zero as we found in this toy example. Instead, we would be satisfied to see that each element of the difference above is comparable in size to machine precision.

``````

(function-lufact-python)=
`````{dropdown} LU factorization (not stable)
:open: true
```{literalinclude} ../python/pkg/FNC/FNC02.py
:filename: lufact.py
:start-line: 31
:end-line: 43
:linenos: true
:language: python
```

```{admonition} About the code
:class: dropdown
Line 5 of {numref}`Function {number} <function-lufact>` points out a subtle issue. Array variables are really just references to blocks of memory. Such a reference is much more efficient to pass around than the complete contents of the array. However, it means that a statement `A_k = A` would just clone the array reference of `A` into the new variable. Any changes made to entries of `A_k` would then also be made to entries of `A`, because they refer to the same location in memory. In this context, we don't want to change the original matrix, so we use `copy` here to create an independent copy of the array contents and a new reference to them.
```
`````
`````

(demo-lu-solve-python)=
``````{dropdown} Solving a linear system by LU factors
:open:

Here are the data for a linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$. 

```{code-cell}
A = np.array([
    [2, 0, 4, 3], 
    [-4, 5, -7, -10], 
    [1, 15, 2, -4.5],
    [-2, 0, 2, -13]
    ])
b = np.array([4, 9, 9, 4])
```

We apply {numref}`Function {number} <function-lufact>` and then do two triangular solves.

```{code-cell}
L, U = FNC.lufact(A)
z = FNC.forwardsub(L, b)
x = FNC.backsub(U, z)
```

A check on the residual assures us that we found the solution.

```{code-cell}
b - A @ x
```

``````


# Example 2.4.3

```{code-cell} ipython3
A = array([[2, 0, 4, 3],[-4, 5, -7, -10],[1, 15, 2, -4.5],[-2, 0, 2, -13]])
```

```{code-cell} ipython3
L,U = FNC.lufact(A)
print("L=",L)
print("U=",U)
```

It's best to compare two floating-point quantities by taking their difference.

```{code-cell} ipython3
LtimesU = L@U
print(LtimesU - A)
```

(Usually we can expect "zero" only up to machine precision. However, all the exact numbers in this example are also floating-point numbers.)

To solve a linear system, we no longer need the matrix $A$. 

```{code-cell} ipython3
b = [4,9,29,40]
z = FNC.forwardsub(L,b)
x = FNC.backsub(U,z)
print(x)
```

```{code-cell} ipython3
print(b - A@x)
```

# Example 2.5.3

```{code-cell} ipython3
n = 6
A = rand(n,n)
x = ones(n)
y = zeros(n)
for i in range(0,n):
    for j in range(0,n):
        y[i] += A[i,j]*x[j]   # 2 flops
```

Each of the loops implies a summation of flops. The total flop count for this algorithm is
\[ \sum_{i=1}^n \sum_{j=1}^n 2 = \sum_{i=1}^n 2n = 2n^2. \]
Since the matrix $A$ has $n^2$ elements, all of which have to be involved in the product, it seems unlikely that we could get a flop count that is smaller than $O(n^2)$.

+++

Let's run an experiment with the built-in matrix-vector multiplication. We assume that flops dominate the computation time and thus measure elapsed time. 

```{code-cell} ipython3
n = 400*arange(1,11)
print("  n           t")
for i,n in enumerate(n):
    A = randn(n,n)  
    x = randn(n)
    start = timer()
    for j in range(1,20): A@x
    print(f"{n:5}   {timer() - start:10.3e}")
```

The reason for doing multiple repetitions at each value of $n$ above is to avoid having times so short that the resolution of the timer is a factor.

+++

# Example 2.5.4

+++

Let's repeat the experiment of the previous example for more, and larger, values of $n$.

```{code-cell} ipython3
N = arange(400,6200,200)
t = zeros(len(N))
for i,n in enumerate(N):
    A = randn(n,n)  
    x = randn(n)
    start = timer()
    for j in range(1,20): A@x
    t[i] = timer() - start
```

Plotting the time as a function of $n$ on log-log scales is equivalent to plotting the logs of the variables, but is formatted more neatly. 

```{code-cell} ipython3
fig,ax = subplots()
ax.loglog(N,t,"-o",label="observed")
ylabel("elapsed time (sec)");
xlabel("$n$");
title("Timing of matrix-vector multiplications");
```

You can see that while the full story is complicated, the graph is trending to a straight line of positive slope. For comparison, we can plot a line that represents $O(n^2)$ growth exactly. (All such lines have slope equal to 2.)

```{code-cell} ipython3
ax.loglog(N,(N/N[-1])**2*t[-1],"--",label="$O(n^2)$")
ax.legend();  fig
```

# Example 2.5.5

+++

We'll test the conclusion of $O(n^3)$ flops experimentally, using the built-in `lu` function instead of the purely instructive `lufact`.

```{code-cell} ipython3
N = arange(200,2600,200)
t = zeros(len(N))
for i,n in enumerate(N):
    A = randn(n,n)  
    start = timer()
    for j in range(1,5): lu(A)
    t[i] = timer() - start
```

We plot the timings on a log-log graph and compare it to $O(n^3)$. The result could vary significantly from machine to machine. 

```{code-cell} ipython3
loglog(N,t,"-o",label="obseved")
loglog(N,(N/N[-1])**3*t[-1],"--",label="$O(n^3)$")
legend();
xlabel("$n$");
ylabel("elapsed time (sec)");
title("Timing of LU factorizations");
```

# Example 2.6.1

+++

Here is the previously solved system.

```{code-cell} ipython3
A = array([[2, 0, 4, 3],[-4, 5, -7, -10],[1, 15, 2, -4.5],[-2, 0, 2, -13]])
b = array([ 4, 9, 29, 40 ])
```

It has a perfectly good solution, obtainable through LU factorization.

```{code-cell} ipython3
L,U = FNC.lufact(A)
x = FNC.backsub( U, FNC.forwardsub(L,b) )
print(x)
```

If we swap the second and fourth equations, nothing essential is changed, and the built-in method still finds the solution.

```{code-cell} ipython3
A[[1,3],:] = A[[3,1],:]  
b[[1,3]] = b[[3,1]]
x = solve(A,b)
print(x)
```

However, LU factorization fails.

```{code-cell} ipython3
L,U = FNC.lufact(A)
print(L)
```

# Example 2.6.2

+++

Here is the system that "broke" LU factorization for us.

```{code-cell} ipython3
A = array([[2, 0, 4, 3],[-4, 5, -7, -10],[1, 15, 2, -4.5],[-2, 0, 2, -13]])
b = array([ 4, 9, 29, 40 ])
```

When we use the `lu` function with three outputs, we get the elements of the PLU factorization.

```{code-cell} ipython3
P,L,U = lu(A)
print(L)
```

```{code-cell} ipython3
print(U)
```

```{code-cell} ipython3
print(P)
```

However, while our notation defines $A=P^TLU$, Python is using $A=PLU$. In order to use the factors to solve the linear system, we have to use $P^T$, which equals $P^{-1}$ for a permutation matrix. 

```{code-cell} ipython3
x = FNC.backsub( U, FNC.forwardsub(L,dot(P.T,b)) )
print(x)
```

If you have to solve many different linear systems for the same matrix, you can perform the computationally expensive factorization just once, and repeat only the much faster triangular solves for the different right-hand sides. 

+++

# Example 2.7.1

+++

The `norm` command computes vector norms.

```{code-cell} ipython3
x = array([2,-3,1,-1])
print( norm(x) )       # 2-norm by default
```

```{code-cell} ipython3
print( norm(x,Inf) )
```

```{code-cell} ipython3
print( norm(x,1) )
```

# Example 2.7.2

```{code-cell} ipython3
A = array([[2,0],[1,-1]])
```

The default matrix norm is *not* the 2-norm. Instead you must provide the 2 explicitly. 

```{code-cell} ipython3
print( norm(A,2) )
```

You can get the 1-norm as well.

```{code-cell} ipython3
print( norm(A,1) )
```

The 1-norm is equivalent to 

```{code-cell} ipython3
max( sum(abs(A),0) )   # sum down the rows (dimension=0)
```

Similarly, we can get the $\infty$-norm and check our formula for it.

```{code-cell} ipython3
print( norm(A,Inf) )
```

```{code-cell} ipython3
max( sum(abs(A),1) )   # sum across columns (dimension=1)
```

Here we illustrate the geometric interpretation of the 2-norm. First, we will sample a lot of vectors on the unit circle in $\mathbb{R}^2$. 

```{code-cell} ipython3
theta = linspace(0,2*pi,601)
x = vstack([cos(theta),sin(theta)])  # 601 unit columns
```

We can apply `A` to every column of `x` simply by using a matrix multiplication.

```{code-cell} ipython3
y = A@x
```

We superimpose the image of the unit circle with the circle whose radius is $\|A\|_2$, and display the plots side by side.

```{code-cell} ipython3
subplot(1,2,1)
plot(x[0,:],x[1,:])
axis("equal")
title("Unit circle")
xlabel("$x_1$")
ylabel("$x_2$");

subplot(1,2,2)
plot(y[0,:],y[1,:])
plot(norm(A,2)*x[0,:],norm(A,2)*x[1,:],"--")
axis("equal")
title("Image under map")
xlabel("$y_1$")
ylabel("$y_2$");
```

# Example 2.8.1

+++

The function `cond` to computes matrix condition numbers. By default, the 2-norm is used. As an example, the family of *Hilbert matrices* is famously badly conditioned. Here is the $7\times 7$  case. 

```{code-cell} ipython3
A = array([[ 1/(i+j+2) for j in range(7)] for i in range(7) ])
print(A)
```

```{code-cell} ipython3
kappa = cond(A)
print(kappa)
```

Next we engineer a linear system problem to which we know the exact answer.

```{code-cell} ipython3
x_exact = 1.0 + arange(7)
b = A@x_exact
```

Now we perturb the data randomly with a vector of norm $10^{-12}$. 

```{code-cell} ipython3
dA = randn(7,7)
dA = 1e-12*(dA/norm(dA,2))
db = randn(7)
db = 1e-12*(db/norm(db,2))
```

We solve the perturbed problem using built-in pivoted LU and see how the solution was changed.

```{code-cell} ipython3
x = solve(A+dA,b+db) 
dx = x - x_exact
```

Here is the relative error in the solution.

```{code-cell} ipython3
print("relative error:",norm(dx) / norm(x_exact))
```

And here are upper bounds predicted using the condition number of the original matrix. 

```{code-cell} ipython3
print("b_bound:",kappa*1e-12/norm(b))
print("A_bound:",kappa*1e-12/norm(A,2))
```

Even if we don't make any manual perturbations to the data, machine epsilon does when we solve the linear system numerically.

```{code-cell} ipython3
x = solve(A,b)
print("relative error:",norm(x-x_exact) / norm(x_exact))
print("rounding bound:",kappa/2**52)
```

Because $\kappa\approx 10^8$, it's possible to lose 8 digits of accuracy in the process of passing from $A$ and $b$ to $x$. That's independent of the algorithm; it's inevitable once the data are expressed in double precision. 

Larger Hilbert matrices are even more poorly conditioned.

```{code-cell} ipython3
A = array([[ 1/(i+j+2) for j in range(14)] for i in range(14) ])
kappa = cond(A)
print(kappa)
```

Before we compute the solution, note that $\kappa$ exceeds `1/eps`. In principle we therefore might end up with an answer that is completely wrong (i.e., a relative error greater than 100%).

```{code-cell} ipython3
print(kappa/2**52)
```

```{code-cell} ipython3
x_exact = 1.0 + arange(14)
b = A@x_exact  
x = solve(A,b)
```

We got an answer. But in fact the error does exceed 100%.

```{code-cell} ipython3
print("relative error:",norm(x-x_exact) / norm(x_exact))
```

# Example 2.9.1

+++

Here is a matrix with both lower and upper bandwidth equal to one. Such a matrix is called tridiagonal.

```{code-cell} ipython3
A = array([ [2, -1,  0,  0,  0,  0],
            [4,  2, -1,  0,  0,  0],
            [0,  3,  0, -1,  0,  0],
            [0,  0,  2,  2, -1,  0],
            [0,  0,  0,  1,  1, -1],
            [0,  0,  0,  0,  0,  2 ]])
```

We can extract the elements on any diagonal using the `diag` command. The "main" or central diagonal is numbered zero, above and to the right of that is positive, and below and to the left is negative.

```{code-cell} ipython3
print( diag(A) )
```

```{code-cell} ipython3
print( diag(A,1) )
```

```{code-cell} ipython3
print( diag(A,-1) )
```

We can also construct matrices by specifying a diagonal with the `diag` function.

```{code-cell} ipython3
A = A + diag([pi,8,6,7],2)
print(A)
```

```{code-cell} ipython3
L,U = FNC.lufact(A)
print(L)
```

```{code-cell} ipython3
print(U)
```

Observe above that the lower and upper bandwidths of $A$ are preserved in the $L$ and $U$ results.

+++

# Example 2.9.2

+++

We'll use a large banded matrix to observe the speedup possible in LU factorization. If we use an ordinary "dense" matrix, though, then there's no way to exploit a banded structure such as tridiagonality.

```{code-cell} ipython3
n = 8000
main = 1 + arange(n)
plusone = linspace(n-1,1,n-1)
minusone = ones(n-1)
A = diag(main) + diag(plusone,1) + diag(minusone,1)
```

```{code-cell} ipython3
start = timer()
lu(A)
print("time:",timer()-start)
```

If instead we construct a proper "sparse" matrix, though, the speedup can be dramatic.

```{code-cell} ipython3
A = sparse.diags([main,plusone,minusone],[0,1,-1],format="csc")
start = timer()
splu(A)
print("time:",timer()-start)
```

# Example 2.9.3

+++

We begin with a symmetric $A$. 

```{code-cell} ipython3
A = array([
    [2,     4,     4,     2],
    [4,     5,     8,    -5],
    [4,     8,     6,     2],
    [2,    -5,     2,   -26]
])
```

Carrying out our usual elimination in the first column leads us to 

```{code-cell} ipython3
L1 = eye(4)
L1[1:,0] = [-2,-2,-1]
A1 = L1@A
print(A1)
```

But now let's note that if we transpose this result, we have the same first column as before! So we could apply  again and then transpose back.

```{code-cell} ipython3
A2 = (L1@A1.T).T
print(A2)
```

Using transpose identities, this is just

```{code-cell} ipython3
A2 = A1@L1.T
print(A2)
```

Now you can see how we proceed down and to the right, eliminating in a column and then symmetrically in the corresponding row.

```{code-cell} ipython3
L2 = eye(4)
L2[2:,1] = [0,-3]
A3 = L2@A2@L2.T
print(A3)
```

Finally, we arrive at a diagonal matrix.

```{code-cell} ipython3
L3 = eye(4)
L3[3,2] = -1
D = L3@A3@L3.T
print(D)
```

# Example 2.9.4

+++

A randomly chosen matrix is extremely unlikely to be symmetric. However, there is a simple way to symmetrize one.

```{code-cell} ipython3
A = 1.0 + floor(9*rand(4,4))
B = A + A.T
print(B)
```

Similarly, a random symmetric matrix is unlikely to be positive definite. The Cholesky algorithm always detects a non-PD matrix by quitting with an error.

```{code-cell} ipython3
cholesky(B)
```

It's not hard to manufacture an SPD matrix to try out the Cholesky factorization.

```{code-cell} ipython3
B = A.T@A
R = cholesky(B)
print(R)
```

```{code-cell} ipython3
norm(R@R.T - B)
```
