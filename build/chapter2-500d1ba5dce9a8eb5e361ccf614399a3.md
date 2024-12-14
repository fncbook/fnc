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

(demo-matrices-julia)=
``````{dropdown} Matrix operations
:open: true

```{note}
While NumPy does have distinct representations for matrices and 2D arrays, use of the explicit matrix class is officially discouraged. We follow this advice here and use arrays to represent both matrices and vectors.
```

:::{index} ! Python; array, ! Python; shape
:::

:::{index}
see: Python; size, Python; shape
:::

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

```{index} 
see: Python; matrix multiplication, Python; \@
```

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

```{index} ! Python; elementwise multiplication; Python; broadcasting
```

The multiplication operator `*` is reserved for elementwise multiplication. Both operands have to be the same size, after any potential broadcasts.

```{code-cell} ipython3
:tags: raises-exception
print(B * A)    # not the same size, so it's an error
```

```{code-cell} ipython3
print((A / 2) * A)
```

```{danger}
If `A` is a matrix, `A**2` is *not* the same as mathematically raising it to the power 2.
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

```{index} Python; broadcasting
```

Most of the mathematical functions, such as cos, sin, log, exp and sqrt, expecting scalars as operands will be broadcast to arrays.

```{code-cell} ipython3
print(np.cos(np.pi * x))      
```
``````


# Example 2.3.2

+++

For a square matrix $A$, the command `solve(A,B)` is mathematically equivalent to $A^{-1}b$. 

```{code-cell} ipython3
A = array([[1,0,-1],[2,2,1],[-1,-3,0]])
b = array([1,2,3])
```

```{code-cell} ipython3
x = solve(A,b)
print(x)
```

One way to check the answer is to compute a quantity known as the **residual**. It is (hopefully) close to machine precision, scaled by the size of the entries of the data. Notice that the matrix–vector multiplication here is performed by `dot`.

```{code-cell} ipython3
residual = b - A@x
print(residual)
```

If the matrix $A$ is singular, you may get an error.

```{code-cell} ipython3
A = array([[0,1],[0,0]])
b = array([1,-1])
solve(A,b)
```

Detecting singularity is a lot like checking whether two floating point numbers are *exactly* equal: because of roundoff, it could be missed. We're headed toward a more robust way to fully describe the situation.

+++

# Example 2.3.3

+++

It's easy to get just the lower triangular part of any matrix using the `tril` command.

```{code-cell} ipython3
A = 1 + floor(9*rand(5,5))
L = tril(A)
print(L)
```

We'll set up and solve a linear system with this matrix.

```{code-cell} ipython3
b = ones(5)
x = FNC.forwardsub(L,b)
print(x)
```

It's not clear what the error in this answer is. However, we can always check the residual.

```{code-cell} ipython3
b - L@x
```

Next we'll engineer a problem to which we know the exact answer. 

```{code-cell} ipython3
alpha = 0.3;
beta = 2.2;
U = diag(ones(5)) + diag([-1,-1,-1,-1],k=1)
U[0,3:5] = [ alpha-beta, beta ]
print(U)
```

```{code-cell} ipython3
x_exact = ones(5)
b = array([alpha,0,0,0,1])

x = FNC.backsub(U,b)
print("error:",x - x_exact)
```

Everything seems OK here. But another example, with a different value for $\beta$, is more troubling.

```{code-cell} ipython3
alpha = 0.3;
beta = 1e12;
U = diag(ones(5)) + diag([-1,-1,-1,-1],k=1)
U[0,[3,4]] = [ alpha-beta, beta ]
b = array([alpha,0,0,0,1])

x = FNC.backsub(U,b)
print("error:",x - x_exact)
```

It's not so good to get four digits of accuracy after starting with sixteen! But the source of the error is not hard to track down. Solving for $x_1$ performs $(\alpha-\beta)+\beta$ in the first row. Since $|\alpha|$ is so much smaller than $|\beta|$, this a recipe for losing digits to subtractive cancellation.

+++

# Example 2.4.1

+++

We create a 4-by-4 linear system with the matrix

```{code-cell} ipython3
A = array([
     [2,    0,    4,    3], 
     [-4,    5,   -7,  -10], 
     [1,   15,    2,   -4.5],
     [-2,    0,    2,  -13]
        ])
```

and with the right-hand side

```{code-cell} ipython3
b = array([ 4, 9, 29, 40 ])
```

We define an augmented matrix by tacking $b$ on the end as a new column.

```{code-cell} ipython3
S = hstack([A,b.reshape(4,1)])
print(S)
```

The goal is to introduce zeros into the lower triangle of this matrix. By using only elementary row operations, we ensure that the matrix $S$ always represents a linear system that is equivalent to the original. We proceed from left to right and top to bottom. The first step is to put a zero in the (2,1) location using a multiple of row 1:

```{code-cell} ipython3
mult21 = S[1,0]/S[0,0]
S[1,:] -= mult21*S[0,:]   # -= means "subtract the following from"
print(S)
```

We repeat the process for the (3,1) and (4,1) entries.

```{code-cell} ipython3
mult31 = S[2,0]/S[0,0];
S[2,:] -= mult31*S[0,:];
mult41 = S[3,0]/S[0,0];
S[3,:] -= mult41*S[0,:];
print(S)
```

The first column has the zero structure we want. To avoid interfering with that, we no longer add multiples of row 1 to anything. Instead, to handle column 2, we use multiples of row 2. We'll also exploit the highly repetitive nature of the operations to write them as a loop. 

```{code-cell} ipython3
for i in range(2,4):
    mult = S[i,1]/S[1,1]
    S[i,:] -= mult*S[1,:]
print(S)
```

We finish out the triangularization with a zero in the (4,3) place. It's a little silly to use a loop for just one iteration, but the point is to establish a pattern.

```{code-cell} ipython3
for i in range(3,4):
    mult = S[i,2]/S[2,2]
    S[i,:] -= mult*S[2,:]
print(S)
```

Recall that $S$ is an augmented matrix: it represents the system $Ux=z$, where

```{code-cell} ipython3
U = S[:,:-1]
z = S[:,-1:].flatten() # turn it into a proper vector
print("U=",U)
print("z=",z)
```

The solutions to this system are identical to those of the original system, but this one can be solved by backward substitution.

```{code-cell} ipython3
x = FNC.backsub(U,z)
print(x)
```

```{code-cell} ipython3
print(b - A@x)
```

# Example 2.4.2

+++

We revisit the previous example using algebra to express the row operations on $A$.

```{code-cell} ipython3
A = array([[2, 0, 4, 3],[-4, 5, -7, -10],[1, 15, 2, -4.5],[-2, 0, 2, -13]])
```

We use the identity and its columns heavily.

```{code-cell} ipython3
I = eye(4)
print(I)
```

The first step is to put a zero in the (2,1) location using a multiple of row 1:

```{code-cell} ipython3
mult21 = A[1,0]/A[0,0]
L21 = I - mult21*outer(I[:,1],I[:,0])
A = dot(L21,A)
print(A)
```

```{code-cell} ipython3
print(L21)
```

We repeat the process for the (3,1) and (4,1) entries. 

```{code-cell} ipython3
mult31 = A[2,0]/A[0,0]
L31 = I - mult31*outer(I[:,2],I[:,0])
A = dot(L31,A)

mult41 = A[3,0]/A[0,0];
L41 = I - mult41*outer(I[:,3],I[:,0])
A = dot(L41,A)

print(A)
```

And so on, following the pattern as before. 

+++

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
