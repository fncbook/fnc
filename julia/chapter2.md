---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---

# Chapter 2 

```{code-cell}
:tags: [remove-cell]
import Pkg; Pkg.activate("/Users/driscoll/Documents/GitHub/fnc")
using FundamentalsNumericalComputation
FNC.init_format()
```

Julia implementations

## Functions

(function-forwardsub-julia)=
``````{dropdown} Forward substitution
:open:
```{literalinclude} ../julia/package/src/chapter02.jl
:filename: forwardsub.jl
:start-line: 0
:end-line: 15
:linenos: true
:language: julia
```

```{index} Julia; sum
```

```{admonition} About the code
:class: dropdown
The `sum` in line 12 gives an error if `i` equals 1, so that case is taken care of before the loop starts.
```
``````

(function-backsub-julia)=
``````{dropdown} Backward substitution
```{literalinclude} ../julia/package/src/chapter02.jl
:filename: backsub.jl
:start-line: 17
:end-line: 32
:linenos: true
:language: julia
```
``````

(function-lufact-julia)=
`````{dropdown} LU factorization (not stable)
```{literalinclude} ../julia/package/src/chapter02.jl
:filename: lufact.jl
:start-line: 34
:end-line: 54
:linenos: true
:language: julia
```

```{admonition} About the code
:class: dropdown
Line 11 of {numref}`Function {number} <function-lufact>` points out two subtle Julia issues. First, vectors and matrix variables are really just references to blocks of memory. Such a reference is much more efficient to pass around than the complete contents of the array. However, it means that a statement `Aₖ=A` just clones the array reference of `A` into the variable `Aₖ`. Any changes made to entries of `Aₖ` would then also be made to entries of `A` because they refer to the same location in memory. In this context we don't want to change the original matrix, so we use `copy` here to create an independent copy of the array contents and a new reference to them.

The second issue is that even when `A` has all integer entries, its LU factors may not. So we convert `Aₖ` to floating point so that line 17 will not fail due to the creation of floating-point values in an integer matrix. An alternative would be to require the caller to provide a floating-point array in the first place.
```
`````

(function-plufact-julia)=
``````{dropdown} LU factorization with partial pivoting
```{literalinclude} ../julia/package/src/chapter02.jl
:filename: plufact.jl
:start-line: 56
:end-line: 80
:linenos: true
:language: julia
```
```
``````

## Examples
### 2.1 @section-linsys-polyinterp
(demo-interp-vander-julia)=
``````{dropdown} @demo-interp-vander
:open: false


We create two vectors for data about the population of China. The first has the years of census data and the other has the population, in millions of people.

```{code-cell}
year = [1982, 2000, 2010, 2015]; 
pop = [1008.18, 1262.64, 1337.82, 1374.62];
```

:::{index} ! Julia; .-, ! Julia; .+
:::

:::{index} Julia; broadcasting
:::

::::{grid} 1 1 2 2
It's convenient to measure time in years since 1980. We use `.-` to subtract a scalar from every element of a vector. We will also use a floating-point value in the subtraction, so the result is also in double precision.
:::{card}
A dotted operator such as `.-` or `.*` acts elementwise, broadcasting scalar values to match up with elements of an array.
:::
::::

```{code-cell}
t = year .- 1980.0
y = pop;
```

:::{index} ! Julia; comprehension
:::

::::{grid} 1 1 2 2
Now we have four data points $(t_1,y_1),\dots,(t_4,y_4)$, so $n=4$ and we seek an interpolating cubic polynomial. We construct the associated Vandermonde matrix:
:::{card}
An expression inside square brackets and ending with a `for` statement is called a **comprehension**. It's often an easy and readable way to construct vectors and matrices. 
:::
::::

```{code-cell}
V = [ t[i]^j for i=1:4, j=0:3 ]
```

:::{index} ! Julia; \\
:::

::::{grid} 1 1 2 2
To solve for the vector of polynomial coefficients, we use a backslash to solve the linear system:
:::{card}
A **backslash** `\` is used to solve a linear system of equations.
:::
::::

```{code-cell}
c = V \ y
```

The algorithms used by the backslash operator are the main topic of this chapter. As a check on the solution, we can compute the *residual*.

```{code-cell} julia
y - V*c
```

Using floating-point arithmetic, it is not realistic to expect exact equality of quantities; a relative difference comparable to $\macheps$ is all we can look for.

::::{grid} 1 1 2 2
By our definitions, the elements of `c` are coefficients in ascending-degree order for the interpolating polynomial. We can use the polynomial to estimate the population of China in 2005:
:::{card}
The `Polynomials` package has functions to make working with polynomials easy and efficient.
:::
::::

```{code-cell}
p = Polynomial(c)    # construct a polynomial
p(2005-1980)         # include the 1980 time shift
```

The official population value for 2005 was 1303.72, so our result is rather good. 

:::{index} ! Julia; scatter
:::

::::{grid} 1 1 2 2
We can visualize the interpolation process. First, we plot the data as points.
:::{card}
The `scatter` function creates a scatter plot of points; you can specify a line connecting the points as well.
::::

```{code-cell}
scatter(t,y, label="actual", legend=:topleft,
    xlabel="years since 1980", ylabel="population (millions)", 
    title="Population of China")
```

:::{index} Julia; range
:::

```{index} ! Julia; broadcasting
```

::::{grid} 1 1 2 2
We want to superimpose a plot of the polynomial. We do that by evaluating it at a vector of points in the interval. The dot after the name of the polynomial is a universal way to apply a function to every element of an array, a technique known as **broadcasting**.
:::{card}
The `range` function constructs evenly spaced values given the endpoints and either the number of values, or the step size between them.

Adding a dot to the end of a function name causes it to be broadcast, i.e., applied to every element of a vector or matrix.
:::
::::

```{code-cell}
# Choose 500 times in the interval [0,35].
tt = range(0,35,length=500)   
# Evaluate the polynomial at all the vector components.
yy = p.(tt)
foreach(println,yy[1:4])
```

:::{index} ! Julia; \!
:::

::::{grid} 1 1 2 2
Now we use `plot!` to add to the current plot, rather than replacing it.
:::{card}
The `plot` function plots lines connecting the given $x$ and $y$ values; you can also specify markers at the points.

By convention, functions whose names end with the bang `!` change the value or state of something, in addition to possibly returning output.
:::
::::

```{code-cell}
plot!(tt,yy,label="interpolant")
```
``````

### 2.2 @section-linsys-matrices
(demo-matrices-julia)=
``````{dropdown} @demo-matrices
:open: false

:::{index} ! Julia; size, ! Julia; length
:::

::::{grid} 1 1 2 2
In Julia, vectors and matrices are one-dimensional and two-dimensional arrays, respectively. Square brackets are used to enclose elements of a matrix or vector. Use spaces for horizontal concatenation, and semicolons or new lines to indicate vertical concatenation.
:::{card}
The `size` function returns the number of rows and columns in a matrix. Use `length` to get the number of elements in a vector or matrix.
:::
::::


```{code-cell}
A = [ 1 2 3 4 5; 50 40 30 20 10
    π sqrt(2) exp(1) (1+sqrt(5))/2 log(3) ]
```

```{code-cell}
m, n = size(A)
```
A vector is not quite the same thing as a matrix: it has only one dimension, not two. Separate its elements by commas or semicolons:

```{code-cell}
x = [ 3, 3, 0, 1, 0 ]
size(x)
```

For some purposes, however, an $n$-vector in Julia is treated like having a column shape. Note the difference if we use spaces instead of commas inside the brackets:

```{code-cell}
y = [ 3 3 0 1 0 ]
size(y)
```

This $1\times 5$ matrix is not equivalent to a vector.

Concatenated elements within brackets may be matrices or vectors for a block representation, as long as all the block sizes are compatible.

```{code-cell}
[ x  x ]
```

```{code-cell}
[ x; x ]
```

The `zeros` and `ones` functions construct matrices with entries all zero or one, respectively.

```{code-cell}
B = [ zeros(3, 2) ones(3, 1) ]
```

```{index} ! Julia; transpose, ! Julia; adjoint, ! Julia; \'
```

A single quote `'` after a matrix returns its adjoint. For real matrices, this is the transpose; for complex-valued matrices, the elements are also conjugated. 

```{code-cell}
A'
```

If `x` is simply a vector, then its transpose has a row shape.

```{code-cell}
x'
```

```{index} ! Julia; range, ! Julia; \:
```

There are many convenient shorthand ways of building vectors and matrices other than entering all of their entries directly or in a loop. To get a range with evenly spaced entries between two endpoints, you have two options. One is to use a colon `:`.

```{code-cell}
y = 1:4              # start:stop
```

```{code-cell}
z = 0:3:12           # start:step:stop
```

(Ranges are not strictly considered vectors, but they behave identically in most circumstances.) Instead of specifying the step size, you can give the number of points in the range if you use `range`.

```{code-cell}
s = range(-1, 1, 5)
```

:::{index} ! Julia; end, ! Julia; indexing arrays
:::

::::{grid} 1 1 2 2
Accessing an element is done by giving one (for a vector) or two (for a matrix) index values within square brackets. 
:::{card}
The `end` keyword refers to the last element in a dimension. It saves you from having to compute and store the size of the matrix first.
::::

```{code-cell}
a = A[2, end-1]
```

```{code-cell}
x[2]
```

The indices can be vectors or ranges, in which case a block of the matrix is accessed.

```{code-cell}
A[1:2, end-2:end]    # first two rows, last three columns
```

```{index} Julia; \:
```

If a dimension has only the index `:` (a colon), then it refers to all the entries in that dimension of the matrix.

```{code-cell}
A[:, 1:2:end]        # all of the odd columns
```

:::{index} ! Julia; diagm
:::

::::{grid} 1 1 2 2
The matrix and vector senses of addition, subtraction, scalar multiplication, multiplication, and power are all handled by the usual symbols. 
:::{card}
Use `diagm` to construct a matrix by its diagonals. A more general syntax puts elements on super- or subdiagonals.
:::
::::

```{code-cell}
B = diagm( [-1, 0, -5] )   # create a diagonal matrix
```

```{code-cell}
@show size(A), size(B);
BA = B * A     # matrix product
```
::::{grid} 1 1 2 2
`A * B` causes an error here, because the dimensions aren't compatible.
:::{card}
Errors are formally called *exceptions* in Julia.
:::
::::

```{code-cell} julia
:tags: raises-exception
A * B    # throws an error
```

A square matrix raised to an integer power is the same as repeated matrix multiplication.

```{code-cell}
B^3    # same as B*B*B
```

Sometimes one instead wants to treat a matrix or vector as a mere array and simply apply a single operation to each element of it. For multiplication, division, and power, the corresponding operators start with a dot.

```{code-cell}
C = -A;
```

Because both matrices are $3\times 5$, `A * C` would be an error here, but elementwise operations are fine.

```{code-cell}
elementwise = A .* C
```

```{index} Julia; broadcasting
```

The two operands of a dot operator have to have the same size—unless one is a scalar, in which case it is expanded or *broadcast* to be the same size as the other operand.

```{code-cell}
x_to_two = x.^2
```

```{code-cell}
two_to_x = 2 .^ x
```

::::{grid} 1 1 2 2
Most of the mathematical functions, such as cos, sin, log, exp, and sqrt, expect scalars as operands. However, you can broadcast any function, including ones that you have defined, across a vector or array by using a special dot syntax.
:::{card}
A dot added to the end of a function name means to apply the function elementwise to an array.
:::
::::

```{code-cell}
show(cos.(π*x))    # broadcast to a function
```

Rather than dotting multiple individual functions, you can use `@.` before an expression to broadcast everything within it.

```{code-cell}
show(@. cos(π*(x+1)^3))    # broadcast an entire expression
```
``````

### 2.3 @section-linsys-linear-systems
(demo-systems-backslash-julia)=
``````{dropdown} @demo-systems-backslash
:open: false
For a square matrix $\mathbf{A}$, the syntax `A \ b` is mathematically equivalent to $\mathbf{A}^{-1} \mathbf{b}$. 

```{code-cell}
A = [1 0 -1; 2 2 1; -1 -3 0]
```

```{code-cell}
b = [1, 2, 3]
```

```{code-cell}
x = A \ b
```

```{index} residual
```

One way to check the answer is to compute a quantity known as the **residual**. It is (ideally) close to machine precision (relative to the elements in the data).

```{code-cell}
residual = b - A*x
```

If the matrix $\mathbf{A}$ is singular, you may get an error.

```{code-cell} julia
:tags: raises-exception
A = [0 1; 0 0]
b = [1, -1]
x = A \ b
```

::::{grid} 1 1 2 2
In this case we can check that the rank of $\mathbf{A}$ is less than its number of columns, indicating singularity.
:::{card}
The function `rank` computes the rank of a matrix. However, it is numerically unstable for matrices that are nearly singular, in a sense to be defined in a later section.
:::
::::

```{code-cell}
rank(A)
```

A linear system with a singular matrix might have no solution or infinitely many solutions, but in either case, backslash will fail. Moreover, detecting singularity is a lot like checking whether two floating-point numbers are *exactly* equal: because of roundoff, it could be missed. In {numref}`section-linsys-condition-number` we'll find a robust way to fully describe this situation.
``````

(demo-systems-triangular-julia)=
``````{dropdown} @demo-systems-triangular
:open: false
```{index} ! Julia; tril, ! Julia; triu
```

::::{grid} 1 1 2 2
It's easy to get just the lower triangular part of any matrix using the `tril` function.
:::{card}
Use `tril` to return a matrix that zeros out everything above the main diagonal. The `triu` function zeros out below the diagonal.
:::
::::

```{code-cell}
A = rand(1.:9., 5, 5)
L = tril(A)
```

We'll set up and solve a linear system with this matrix.

```{code-cell}
b = ones(5)
x = FNC.forwardsub(L,b)
```

It's not clear how accurate this answer is. However, the residual should be zero or comparable to $\macheps$.

```{code-cell}
b - L * x
```

```{index} ! Julia; Pair, Julia; diagm
```

::::{grid} 1 1 2 2
Next we'll engineer a problem to which we know the exact answer. Use `\alpha` <kbd>Tab</kbd> and `\beta` <kbd>Tab</kbd> to get the Greek letters.
:::{card}
The notation `0=>ones(5)` creates a `Pair`. In `diagm`, pairs indicate the position of a diagonal and the elements that are to be placed on it.
:::
::::

```{code-cell}
α = 0.3;
β = 2.2;
U = diagm( 0=>ones(5), 1=>[-1, -1, -1, -1] )
U[1, [4, 5]] = [ α - β, β ]
U
```

```{code-cell}
x_exact = ones(5)
b = [α, 0, 0, 0, 1]
```

Now we use backward substitution to solve for $\mathbf{x}$, and compare to the exact solution we know already.

```{code-cell}
x = FNC.backsub(U,b)
err = x - x_exact
```

Everything seems OK here. But another example, with a different value for $\beta$, is more troubling.

```{code-cell}
α = 0.3;
β = 1e12;
U = diagm( 0=>ones(5), 1=>[-1, -1, -1, -1] )
U[1, [4, 5]] = [ α - β, β ]
b = [α, 0, 0, 0, 1]

x = FNC.backsub(U,b)
err = x - x_exact
```

It's not so good to get 4 digits of accuracy after starting with 16! The source of the error is not hard to track down. Solving for $x_1$ performs $(\alpha-\beta)+\beta$ in the first row. Since $|\alpha|$ is so much smaller than $|\beta|$, this a recipe for losing digits to subtractive cancellation.
``````
### 2.4 @section-linsys-lu
(demo-lu-outertri-julia)= 
``````{dropdown} @demo-lu-outertri
:open: false

```{index} Julia; tril, Julia; triu
```
We explore the outer product formula for two random triangular matrices.

```{code-cell}
L = tril( rand(1:9, 3, 3) )
```

```{code-cell}
U = triu( rand(1:9, 3, 3) )
```

::::{grid} 1 1 2 2
Here are the three outer products in the sum in {eq}`matrixouter`:
:::{card}
Although `U[1,:]` is a row of `U`, it is a vector, and as such it has a default column interpretation.
:::
::::

```{code-cell}
L[:, 1] * U[1, :]'
```

```{code-cell}
L[:, 2] * U[2, :]'
```

```{code-cell}
L[:, 3] * U[3, :]'
```

Simply because of the triangular zero structures, only the first outer product contributes to the first row and first column of the entire product. 
``````

(demo-lu-derive-julia)=
``````{dropdown} @demo-lu-derive
For illustration, we work on a $4 \times 4$ matrix. We name it with a subscript in preparation for what comes.

```{code-cell}
A₁ = [
     2    0    4     3 
    -4    5   -7   -10 
     1   15    2   -4.5
    -2    0    2   -13
    ];
L = diagm(ones(4))
U = zeros(4, 4);
```

Now we appeal to {eq}`outer-row1`. Since $L_{11}=1$, we see that the first row of $\mathbf{U}$ is just the first row of $\mathbf{A}_1$.

```{code-cell}
U[1, :] = A₁[1, :]
U
```

From {eq}`outer-col1`, we see that we can find the first column of $\mathbf{L}$ from the first column of $\mathbf{A}_1$. 

```{code-cell}
L[:, 1] = A₁[:, 1] / U[1, 1]
L
```

 We have obtained the first term in the sum {eq}`matrixouter` for $\mathbf{L}\mathbf{U}$, and we subtract it away from $\mathbf{A}_1$.

```{code-cell}
A₂ = A₁ - L[:, 1] * U[1, :]'
```

Now $\mathbf{A}_2 = \boldsymbol{\ell}_2\mathbf{u}_2^T + \boldsymbol{\ell}_3\mathbf{u}_3^T + \boldsymbol{\ell}_4\mathbf{u}_4^T.$ If we ignore the first row and first column of the matrices in this equation, then in what remains we are in the same situation as at the start. Specifically, only $\boldsymbol{\ell}_2\mathbf{u}_2^T$ has any effect on the second row and column, so we can deduce them now.

```{code-cell}
U[2, :] = A₂[2, :]
L[:, 2] = A₂[:, 2] / U[2, 2]
L
```

If we subtract off the latest outer product, we have a matrix that is zero in the first *two* rows and columns. 

```{code-cell}
A₃ = A₂ - L[:, 2] * U[2, :]'
```

Now we can deal with the lower right $2\times 2$ submatrix of the remainder in a similar fashion.

```{code-cell}
U[3, :] = A₃[3, :]
L[:, 3] = A₃[:, 3] / U[3, 3]
A₄ = A₃ - L[:, 3] * U[3, :]'
```

Finally, we pick up the last unknown in the factors.

```{code-cell}
U[4, 4] = A₄[4, 4];
```

We now have all of $\mathbf{L}$,

```{code-cell} 
L
```

and all of $\mathbf{U}$,

```{code-cell}
U
```

We can verify that we have a correct factorization of the original matrix by computing the backward error:

```{code-cell} 
A₁ - L * U
```

IIn floating point, we cannot expect the difference to be exactly zero as we found in this toy example. Instead, we would be satisfied to see that each element of the difference above is comparable in size to machine precision.

``````

(demo-lu-solve-julia)=
``````{dropdown} @demo-lu-solve
Here are the data for a linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$. 

```{code-cell}
A = [2 0 4 3; -4 5 -7 -10; 1 15 2 -4.5; -2 0 2 -13];
b = [4,9,9,4];
```

We apply {numref}`Function {number} <function-lufact>` and then do two triangular solves.

```{code-cell}
L, U = FNC.lufact(A)
z = FNC.forwardsub(L, b)
x = FNC.backsub(U, z)
```

A check on the residual assures us that we found the solution.

```{code-cell}
b - A*x
```
``````

### 2.5 @section-linsys-efficiency

(demo-flops-mvmult-julia)=
``````{dropdown} @demo-flops-mvmult
Here is a straightforward implementation of matrix-vector multiplication.

```{code-cell}
n = 6
A = randn(n, n)
x = rand(n)
y = zeros(n)
for i in 1:n
    for j in 1:n
        y[i] += A[i, j] * x[j]    # 1 multiply, 1 add
    end
end
```

Each of the loops implies a summation of flops. The total flop count for this algorithm is

$$ \sum_{i=1}^n \sum_{j=1}^n 2 = \sum_{i=1}^n 2n = 2n^2. $$

Since the matrix $\mathbf{A}$ has $n^2$ elements, all of which have to be involved in the product, it seems unlikely that we could get a flop count that is smaller than $O(n^2)$ in general.

```{index} ! Julia; push\!, ! Julia; for
```

::::{grid} 1 1 2 2
Let's run an experiment with the built-in matrix-vector multiplication. Note that Julia is unusual in that loops have a variable scope separate from its enclosing code. Thus, `for n in n` below means that inside the loop, the name `n` will take on each one of the values that were previously assigned to the vector `n`.
:::{card}
The `push!` function attaches a new value to the end of a vector.
:::
::::

```{code-cell}
n = 1000:1000:5000
t = []
for n in n
    A = randn(n, n)  
    x = randn(n)
    time = @elapsed for j in 1:80; A * x; end
    push!(t, time)
end
```

The reason for doing multiple repetitions at each value of $n$ in the loop above is to avoid having times so short that the resolution of the timer is significant.

```{code-cell}
pretty_table([n t], header=(["size", "time"], ["", "(sec)"]))
```

```{index} Julia; Boolean indexing
```

::::{grid} 1 1 2 2
Looking at the timings just for $n=2000$ and $n=4000$, they have ratio
:::{card}
The expression `n.==4000` here produces a vector of Boolean (true/false) values the same size as `n`. This result is used to index within `t`, accessing only the value for which the comparison is true.
:::
::::

```{code-cell}
@show t[n.==4000] ./ t[n.==2000];
```

If the run time is dominated by flops, then we expect this ratio to be

$$
\frac{2(4000)^2}{2(2000)^2}=4.
$$
``````

(demo-flops-loglog-julia)=
``````{dropdown} @demo-flops-loglog
Let's repeat the experiment of the previous figure for more, and larger, values of $n$.

```{code-cell}
randn(5,5)*randn(5);  # throwaway to force compilation

n = 400:200:6000
t = []
for n in n
    A = randn(n, n)  
    x = randn(n)
    time = @elapsed for j in 1:50; A * x; end
    push!(t, time)
end
```

Plotting the time as a function of $n$ on log-log scales is equivalent to plotting the logs of the variables.

```{code-cell}
scatter(n, t, label="data", legend=false,
    xaxis=(:log10, L"n"), yaxis=(:log10, "elapsed time (sec)"),
    title="Timing of matrix-vector multiplications")
```

You can see that while the full story is complicated, the graph is trending to a straight line of positive slope. For comparison, we can plot a line that represents $O(n^2)$ growth exactly. (All such lines have slope equal to 2.)

```{code-cell}
plot!(n, t[end] * (n/n[end]).^2, l=:dash,
    label=L"O(n^2)", legend=:topleft)
```
``````

(demo-flops-lufact-julia)=
``````{dropdown} @demo-flops-lufact
::::{grid} 1 1 2 2
We'll test the conclusion of $O(n^3)$ flops experimentally, using the built-in `lu` function instead of the purely instructive `lufact`.
:::{card}
The first time a function is invoked, there may be significant time needed to compile it in memory. Thus, when timing a function, run it at least once before beginning the timing.
:::
::::

```{code-cell}
lu(randn(3, 3));   # throwaway to force compilation

n = 400:400:4000
t = []
for n in n
    A = randn(n, n)  
    time = @elapsed for j in 1:12; lu(A); end
    push!(t, time)
end
```

We plot the timings on a log-log graph and compare it to $O(n^3)$. The result could vary significantly from machine to machine, but in theory the data should start to parallel the line as $n\to\infty$.

```{code-cell}
scatter(n, t, label="data", legend=:topleft,
    xaxis=(:log10, L"n"), yaxis=(:log10, "elapsed time"))
plot!(n, t[end ]* (n/n[end]).^3, l=:dash, label=L"O(n^3)")
```
``````

### 2.6 @section-linsys-pivoting
(demo-pivoting-fail-julia)=
``````{dropdown} @demo-pivoting-fail
Here is a previously encountered matrix that factors well.

```{code-cell}
A = [2 0 4 3 ; -4 5 -7 -10 ; 1 15 2 -4.5 ; -2 0 2 -13];
L, U = FNC.lufact(A)
L
```

If we swap the second and fourth rows of $\mathbf{A}$, the result is still nonsingular. However, the factorization now fails.

```{code-cell}
A[[2, 4], :] = A[[4, 2], :]  
L, U = FNC.lufact(A)
L
```

```{index} Julia; NaN
```

The presence of `NaN` in the result indicates that some impossible operation was required. The source of the problem is easy to locate. We can find the first outer product in the factorization just fine:

```{code-cell}
U[1, :] = A[1, :]
L[:, 1] = A[:, 1] / U[1, 1]
A -= L[:, 1] * U[1, :]'
```

The next step is `U[2, :] = A[2, :]`, which is also OK. But then we are supposed to divide by `U[2, 2]`, which is zero. The algorithm cannot continue.
``````

(demo-pivoting-fix-julia)=
``````{dropdown} @demo-pivoting-fix
Here is the trouble-making matrix from {numref}`Demo {number} <demo-pivoting-fail>`.

```{code-cell}
A₁ = [2 0 4 3 ; -2 0 2 -13; 1 15 2 -4.5 ; -4 5 -7 -10]
```

::::{grid} 1 1 2 2
We now find the largest candidate pivot in the first column. We don't care about sign, so we take absolute values before finding the max.
:::{card}
The `argmax` function returns the location of the largest element of a vector or matrix.
:::
::::


```{code-cell}
i = argmax( abs.(A₁[:, 1]) ) 
```

This is the row of the matrix that we extract to put into $\mathbf{U}$. That guarantees that the division used to find $\boldsymbol{\ell}_1$ will be valid.

```{code-cell}
L, U = zeros(4,4),zeros(4,4)
U[1, :] = A₁[i, :]
L[:, 1] = A₁[:, 1] / U[1, 1]
A₂ = A₁ - L[:, 1] * U[1, :]'
```

Observe that $\mathbf{A}_2$ has a new zero row and zero column, but the zero row is the fourth rather than the first. However, we forge on by using the largest possible pivot in column 2 for the next outer product.

```{code-cell}
@show i = argmax( abs.(A₂[:, 2]) ) 
U[2, :] = A₂[i, :]
L[:, 2] = A₂[:, 2] / U[2, 2]
A₃ = A₂ - L[:, 2] * U[2, :]'
```

Now we have zeroed out the third row as well as the second column. We can finish out the procedure.

```{code-cell}
@show i = argmax( abs.(A₃[:, 3]) ) 
U[3, :] = A₃[i, :]
L[:, 3] = A₃[:, 3] / U[3, 3]
A₄ = A₃ - L[:, 3] * U[3, :]'
```

```{code-cell}
@show i = argmax( abs.(A₄[:, 4]) ) 
U[4, :] = A₄[i, :]
L[:, 4] = A₄[:, 4] / U[4, 4];
```

We do have a factorization of the original matrix:

```{code-cell}
A₁ - L * U
```

And $\mathbf{U}$ has the required structure:

```{code-cell}
U
```

However, the triangularity of $\mathbf{L}$ has been broken.

```{code-cell}
L
```
``````

(demo-pivoting-permute-julia)=
``````{dropdown} @demo-pivoting-permute
Here again is the matrix from {numref}`Demo {number} <demo-pivoting-fix>`.

```{code-cell}
A = [2 0 4 3 ; -2 0 2 -13; 1 15 2 -4.5 ; -4 5 -7 -10]
```

As the factorization proceeded, the pivots were selected from rows 4, 3, 2, and finally 1. If we were to put the rows of $\mathbf{A}$ into that order, then the algorithm would run exactly like the plain LU factorization from {numref}`section-linsys-lu`. 

```{code-cell}
B = A[[4, 3, 2, 1], :]
L, U = FNC.lufact(B);
```

We obtain the same $\mathbf{U}$ as before:

```{code-cell}
U
```

And $\mathbf{L}$ has the same rows as before, but arranged into triangular order:

```{code-cell}
L
```
``````

(demo-pivoting-usage-julia)=
``````{dropdown} @demo-pivoting-usage
The third output of `plufact` is the permutation vector we need to apply to $\mathbf{A}$.

```{code-cell}
A = rand(1:20, 4, 4)
L, U, p = FNC.plufact(A)
A[p,:] - L * U   # should be ≈ 0
```

Given a vector $\mathbf{b}$, we solve $\mathbf{A}\mathbf{x}=\mathbf{b}$ by first permuting the entries of $\mathbf{b}$ and then proceeding as before.

```{code-cell}
b = rand(4)
z = FNC.forwardsub(L,b[p])
x = FNC.backsub(U,z)
```

A residual check is successful:

```{code-cell}
b - A*x
```
``````

(demo-pivoting-builtin-julia)=
``````{dropdown} @demo-pivoting-builtin
With the syntax `A \ b`, the matrix `A` is PLU-factored, followed by two triangular solves.

```{code-cell}
A = randn(500, 500)   # 500x500 with normal random entries
A \ rand(500)          # force compilation
@elapsed for k=1:50; A \ rand(500); end
```

In {numref}`section-linsys-efficiency` we showed that the factorization is by far the most costly part of the solution process. A factorization object allows us to do that costly step only once per unique matrix. 

```{code-cell}
factored = lu(A)     # store factorization result
factored \ rand(500)   # force compilation
@elapsed for k=1:50; factored \ rand(500); end
```
``````

(demo-pivoting-stable-julia)=
``````{dropdown} @demo-pivoting-stable
We construct a linear system for this matrix with $\epsilon=10^{-12}$ and exact solution $[1,1]$:

```{code-cell}
ϵ = 1e-12
A = [-ϵ 1; 1 -1]
b = A * [1, 1]
```

We can factor the matrix without pivoting and solve for $\mathbf{x}$.

```{code-cell}
L, U = FNC.lufact(A)
x = FNC.backsub( U, FNC.forwardsub(L, b) )
```

Note that we have obtained only about 5 accurate digits for $x_1$. We could make the result even more inaccurate by making $\epsilon$ even smaller:

```{code-cell}
ϵ = 1e-20; A = [-ϵ 1; 1 -1]
b = A * [1, 1]
L, U = FNC.lufact(A)
x = FNC.backsub( U, FNC.forwardsub(L, b) )
```

This effect is not due to ill conditioning of the problem—a solution with PLU factorization works perfectly:

```{code-cell}
A \ b
```
``````

### 2.7 @section-linsys-norms
(demo-norms-vector-julia)=
``````{dropdown} @demo-norms-vector

```{index} ! Julia; norm
```

In Julia the `LinearAlgebra` package has a `norm` function for vector norms.

```{code-cell}
x = [2, -3, 1, -1]
twonorm = norm(x)         # or norm(x,2)
```

```{code-cell}
infnorm = norm(x, Inf)
```

```{code-cell}
onenorm = norm(x, 1)
```

```{index} ! Julia; normalize
```

There is also a `normalize` function that divides a vector by its norm, making it a unit vector.

```{code-cell}
normalize(x, Inf)
```
``````

(demo-norms-matrix-julia)=
``````{dropdown} @demo-norms-matrix
```{code-cell}
A = [ 2 0; 1 -1 ]
```

In Julia, one uses `norm` for vector norms and for the Frobenius norm of a matrix, which is like stacking the matrix into a single vector before taking the 2-norm. 

```{code-cell}
Fronorm = norm(A)
```

```{index} ! Julia; opnorm
```

Most of the time we want to use `opnorm`, which is an induced matrix norm. The default is the 2-norm.

```{code-cell}
twonorm = opnorm(A)
```

You can get the 1-norm as well.

```{code-cell}
onenorm = opnorm(A, 1)
```

```{index} ! Julia; maximum, ! Julia; minimum, ! Julia; sum
```

::::{grid} 1 1 2 2
According to {eq}`mxonenorm`, the matrix 1-norm is equivalent to the maximum of the sums down the columns (in absolute value).
:::{card}
Use `sum` to sum along a dimension of a matrix. You can also sum over the entire matrix by omitting the `dims` argument.

The `maximum` and `minimum` functions also work along one dimension or over an entire matrix. To get both values at once, use `extrema`.
:::
::::

```{code-cell}
# Sum down the rows (1st matrix dimension):
maximum( sum(abs.(A), dims=1) )   
```

Similarly, we can get the $\infty$-norm and check our formula for it.

```{code-cell}
infnorm = opnorm(A, Inf)
```

```{code-cell}
 # Sum across columns (2nd matrix dimension):
maximum( sum(abs.(A), dims=2) )  
```
::::{grid} 1 1 2 2
Next we illustrate a geometric interpretation of the 2-norm. First, we will sample a lot of vectors on the unit circle in $\mathbb{R}^2$.
:::{card}
You can use functions as values, e.g., as elements of a vector. 
:::
::::

```{code-cell}
# Construct 601 unit column vectors.
θ = 2π * (0:1/600:1)   # type \theta then Tab
x = [ fun(t) for fun in [cos, sin], t in θ ];
```

```{index} ! Julia; subplots
```

To create an array of plots, start with a `plot` that has a `layout` argument, then do subsequent `plot!` calls with a `subplot` argument.

```{code-cell}
plot(aspect_ratio=1, layout=(1, 2),
    xlabel=L"x_1",  ylabel=L"x_2")
plot!(x[1, :], x[2, :], subplot=1, title="Unit circle") 
```

The linear function $\mathbf{f}(\mathbf{x}) = \mathbf{A}\mathbf{x}$ defines a mapping from $\mathbb{R}^2$ to $\mathbb{R}^2$. We can apply `A` to every column of `x` by using a single matrix multiplication.

```{code-cell}
Ax = A * x;
```

The image of the transformed vectors is an ellipse. 

```{code-cell}
plot!(Ax[1, :], Ax[2, :], 
    subplot=2, title="Image under x → Ax")
```

That ellipse just touches the circle of radius $\|\mathbf{A}\|_2$.

```{code-cell}
plot!(twonorm*x[1, :], twonorm*x[2, :], subplot=2, l=:dash)
```
``````


### 2.8 @section-linsys-condition-number
(demo-condition-bound-julia)=
``````{dropdown} @demo-condition-bound

```{index} ! Julia; cond
```
::::{grid} 1 1 2 2
Julia has a function `cond` to compute matrix condition numbers. By default, the 2-norm is used. As an example, the family of *Hilbert matrices* is famously badly conditioned. Here is the $6\times 6$  case.
:::{card}
Type `\kappa` followed by <kbd>Tab</kbd> to get the Greek letter $\kappa$.
:::
::::

```{code-cell}
A = [ 1 / (i + j) for i in 1:6, j in 1:6 ]
κ = cond(A)
```

Because $\kappa\approx 10^8$, it's possible to lose nearly 8 digits of accuracy in the process of passing from $\mathbf{A}$ and $\mathbf{b}$ to $\mathbf{x}$. That fact is independent of the algorithm; it's inevitable once the data are expressed in finite precision. 

Let's engineer a linear system problem to observe the effect of a perturbation. We will make sure we know the exact answer.

```{code-cell}
x = 1:6
b = A * x
```

Now we perturb the system matrix and vector randomly by $10^{-10}$ in norm.

```{code-cell} 
# type \Delta then Tab to get Δ
ΔA = randn(size(A));  ΔA = 1e-10 * (ΔA / opnorm(ΔA));
Δb = randn(size(b));  Δb = 1e-10 * normalize(Δb);
```

We solve the perturbed problem using pivoted LU and see how the solution was changed.

```{code-cell}
new_x = ((A + ΔA) \ (b + Δb))
Δx = new_x - x
```

Here is the relative error in the solution.

```{code-cell}
@show relative_error = norm(Δx) / norm(x);
```

And here are upper bounds predicted using the condition number of the original matrix.

```{code-cell}
println("Upper bound due to b: $(κ * norm(Δb) / norm(b))")
println("Upper bound due to A: $(κ * opnorm(ΔA) / opnorm(A))")
```

Even if we didn't make any manual perturbations to the data, machine roundoff does so at the relative level of $\macheps$.

```{code-cell}
Δx = A\b - x
@show relative_error = norm(Δx) / norm(x);
@show rounding_bound = κ * eps();
```

Larger Hilbert matrices are even more poorly conditioned:

```{code-cell}
A = [ 1 / (i + j) for i=1:14, j=1:14 ];
κ = cond(A)
```

Note that $\kappa$ exceeds $1/\macheps$. In principle we therefore may end up with an answer that has relative error greater than 100%.

```{code-cell}
rounding_bound = κ*eps()
```

Let's put that prediction to the test.

```{code-cell}
x = 1:14
b = A * x  
Δx = A\b - x
@show relative_error = norm(Δx) / norm(x);
```

As anticipated, the solution has zero accurate digits in the 2-norm.
``````

### 2.9 @section-linsys-structure
(demo-structure-banded-julia)=
``````{dropdown} @demo-structure-banded
```{index} ! Julia; fill, Julia; diagm, ! Julia; diag
```

::::{grid} 1 1 2 2
Here is a small tridiagonal matrix. Note that there is one fewer element on the super- and subdiagonals than on the main diagonal.
:::{card}
Use `fill` to create an array of a given size, with each element equal to a provided value.
:::
::::

```{code-cell}
A = diagm( -1 => [4, 3, 2, 1, 0], 
    0 => [2, 2, 0, 2, 1, 2], 
    1 => fill(-1, 5) )
```

```{index} ! Julia; diag
```

::::{grid} 1 1 2 2
We can extract the elements on any diagonal using the `diag` function. The main or central diagonal is numbered zero, above and to the right of that is positive, and below and to the left is negative.
:::{card}
The `diag` function extracts the elements from a specified diagonal of a matrix.
:::
::::

```{code-cell}
@show diag_main = diag(A);
@show diag_minusone = diag(A, -1);
```
The lower and upper bandwidths of $\mathbf{A}$ are repeated in the factors from the unpivoted LU factorization. 

```{code-cell}
L, U = FNC.lufact(A)
L
```

```{code-cell}
U
```
``````

(demo-structure-symm-julia)=
``````{dropdown} @demo-structure-symm

We begin with a symmetric $\mathbf{A}$.

```{code-cell}
A₁ = [  2     4     4     2
        4     5     8    -5
        4     8     6     2
        2    -5     2   -26 ];
```

We won't use pivoting, so the pivot element is at position (1,1). This will become the first element on the diagonal of $\mathbf{D}$. Then we divide by that pivot to get the first column of $\mathbf{L}$.

```{code-cell}
L = diagm(ones(4))
d = zeros(4)
d[1] = A₁[1, 1]
L[:, 1] = A₁[:, 1] / d[1]
A₂ = A₁ - d[1] * L[:, 1] * L[:, 1]'
```

We are now set up the same way for the submatrix in rows and columns 2–4.

```{code-cell}
d[2] = A₂[2, 2]
L[:, 2] = A₂[:, 2] / d[2]
A₃ = A₂ - d[2] * L[:, 2] * L[:, 2]'
```

We continue working our way down the diagonal.

```{code-cell}
d[3] = A₃[3, 3]
L[:, 3] = A₃[:, 3] / d[3]
A₄ = A₃ - d[3] * L[:, 3] * L[:, 3]'
d[4] = A₄[4, 4]
@show d;
L
```

We have arrived at the desired factorization, which we can validate:

```{code-cell}
opnorm(A₁ - (L * diagm(d) * L'))
```
``````

(demo-structure-cholesky-julia)=
``````{dropdown} @demo-structure-cholesky
A randomly chosen matrix is extremely unlikely to be symmetric. However, there is a simple way to symmetrize one.

```{code-cell}
A = rand(1.0:9.0, 4, 4)
B = A + A'
```

```{index} ! Julia; cholesky
```

::::{grid} 1 1 2 2
Similarly, a random symmetric matrix is unlikely to be positive definite. The Cholesky algorithm always detects a non-PD matrix by quitting with an error.
:::{card}
The `cholesky` function computes a Cholesky factorization if possible, or throws an error for a non-positive-definite matrix. However, it does *not* check for symmetry.
:::
::::

```{code-cell}
:tags: raises-exception
cholesky(B)    # throws an error
```

It's not hard to manufacture an SPD matrix to try out the Cholesky factorization.

```{code-cell}
B = A' * A
cf = cholesky(B)
```

What's returned is a factorization object. Another step is required to extract the factor as a matrix:

```{code-cell}
R = cf.U
```

Here we validate the factorization:

```{code-cell}
opnorm(R' * R - B) / opnorm(B)
```
``````
