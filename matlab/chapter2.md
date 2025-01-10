---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---

# Chapter 2

MATLAB implementations

## Functions

(function-forwardsub-matlab)=
``````{dropdown} Forward substitution
:open:
```{literalinclude} FNC-matlab/forwardsub.m
:linenos: true
:language: matlab
```

```{admonition} About the code
:class: dropdown
Line 12 implements @forwardsub. It contains an inner product between row $i$ of $\mathbf{L}$ and the solution vector $\mathbf{x}$, using only the entries of $\mathbf{x}$ that have already been computed. 
```
``````

(function-backsub-matlab)=
``````{dropdown} Backward substitution
:open:
```{literalinclude} FNC-matlab/backsub.m
:linenos: true
:language: matlab
```
``````

(function-lufact-matlab)=
`````{dropdown} LU factorization (not stable)
:open:
```{literalinclude} FNC-matlab/lufact.m
:linenos: true
:language: matlab
```
`````

(function-plufact-matlab)=
``````{dropdown} LU factorization with partial pivoting
:open:
```{literalinclude} FNC-matlab/plufact.m
:linenos: true
:language: matlab
```
```
``````


## Examples

```{code-cell}
:tags: [remove-cell]
cd  /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init
pwd
```

### 2.1 @section-linsys-polyinterp
(demo-interp-vander-matlab)=
``````{dropdown} @demo-interp-vander
We create two column vectors for data about the population of China. The first has the years of census data and the other has the population, in millions of people.

```{code-cell}
year = [1982; 2000; 2010; 2015]; 
pop = [1008.18; 1262.64; 1337.82; 1374.62];
```

It's convenient to measure time in years since 1980. 

```{code-cell}
t = year - 1980;
y = pop;
```

```{index} ! MATLAB; vander
```

Now we have four data points $(t_1,y_1),\dots,(t_4,y_4)$, so $n=4$ and we seek an interpolating cubic polynomial. We construct the associated Vandermonde matrix:

```{code-cell}
V = vander(t)
```

:::{index} ! MATLAB; \\
:::

To solve for the vector of polynomial coefficients, we use a backslash to solve the linear system:
```{tip}
:class: dropdown
A **backslash** `\` is used to solve a linear system of equations.
```

```{code-cell}
c = V \ y
```

The algorithms used by the backslash operator are the main topic of this chapter. As a check on the solution, we can compute the *residual*.

```{code-cell}
y - V * c
```

Using floating-point arithmetic, it is not realistic to expect exact equality of quantities; a relative difference comparable to $\macheps$ is all we can look for.

By our definitions, the elements of `c` are coefficients in descending-degree order for the interpolating polynomial. We can use the polynomial to estimate the population of China in 2005:

```{code-cell}
p = @(t) polyval(c, t - 1980);  % include the 1980 time shift
p(2005)
```

The official population value for 2005 was 1303.72, so our result is rather good. 

:::{index} ! MATLAB; scatter
:::

We can visualize the interpolation process. First, we plot the data as points.
```{tip}
:class: dropdown
The `scatter` function creates a scatter plot of points; you can specify a line connecting the points as well.
```

```{code-cell}
scatter(year, y)
xlabel("years since 1980")
ylabel("population (millions)")
title(("Population of China"));
```

:::{index} ! MATLAB; linspace
:::

We want to superimpose a plot of the polynomial. We do that by evaluating it at a vector of points in the interval.
```{tip}
:class: dropdown
The `linspace` function constructs evenly spaced values given the endpoints and the number of values.
```

```{code-cell}
tt = linspace(1980, 2015, 500);    % 500 times in the interval [1980, 2015]
yy = p(tt);                        % evaluate p at all the vector elements
yy(1:4)
```

```{index} ! MATLAB; hold on, ! MATLAB; plot
```

Now we use `plot!` to add to the current plot, rather than replacing it.
```{tip}
:class: dropdown
Use `hold on` to add to an existing plot rather than replacing it.
The `plot` function plots lines connecting the given $x$ and $y$ values; you can also specify markers at the points.
```

```{code-cell}
hold on 
plot(tt, yy)
legend("data", "interpolant", "location", "northwest");
```
``````

### 2.2 @section-linsys-matrices
(demo-matrices-matlab)=
``````{dropdown} @demo-matrices
:::{index} ! MATLAB; size, ! MATLAB; length
:::

In MATLAB, every numerical value is treated like a matrix. A matrix with one row or one column is interpreted as a vector, and a $1\times 1$ matrix is interpreted as a scalar. 

Square brackets are used to enclose elements of a matrix or vector. Use spaces for horizontal concatenation, and semicolons or new lines to indicate vertical concatenation.
```{tip}
:class: dropdown
The `size` function returns the number of rows and columns in a matrix. Use `length` to get the number of elements in a vector or matrix.
```


```{code-cell}
A = [ 
    1       2      3             4      5; 
    50     40     30            20     10
    pi sqrt(2) exp(1) (1+sqrt(5))/2 log(3) 
    ]
```

```{code-cell}
m, n = size(A)
```

```{code-cell}
x = [ 3, 3, 0, 1, 0 ];   % row vector
size(x)
```

Concatenated elements within brackets may be matrices or vectors for a block representation, as long as all the block sizes are compatible.

```{code-cell}
[ x  x ]
```

```{code-cell}
[ x; x ]
```

```{index} ! MATLAB; zeros, ! MATLAB; ones
```

The `zeros` and `ones` functions construct matrices with entries all zero or one, respectively.

```{code-cell}
B = [ zeros(3, 2) ones(3, 1) ]
```

```{index} ! MATLAB; transpose, ! MATLAB; adjoint, ! MATLAB; \'
```

A single quote `'` after a matrix returns its adjoint. For real matrices, this is the transpose; for complex-valued matrices, the elements are also conjugated. 

```{code-cell}
A'
```

```{index} ! MATLAB; linspace, ! MATLAB; \:
```

There are many convenient shorthand ways of building vectors and matrices other than entering all of their entries directly or in a loop. To get a range with evenly spaced entries between two endpoints, you have two options. One is to use a colon `:`.

```{code-cell}
y = 1:4              % start:stop
```

```{code-cell}
z = 0:3:12           % start:step:stop
```

Instead of specifying the step size, you can give the number of points in the range if you use `linspace`.

```{code-cell}
s = linspace(-1, 1, 5)    % row result
```

:::{index} ! MATLAB; end, ! MATLAB; indexing arrays
:::

Accessing an element is done by giving one (for a vector) or two (for a matrix) index values within parentheses. 
```{tip}
:class: dropdown
The `end` keyword refers to the last element in a dimension. It saves you from having to compute and store the size of the matrix first.
```

```{code-cell}
a = A(2, end-1)
```

```{code-cell}
x(2)
```

The indices can be vectors or ranges, in which case a block of the matrix is accessed.

```{code-cell}
A(1:2, end-2:end)    % first two rows, last three columns
```

```{index} MATLAB; \:
```

If a dimension has only the index `:` (a colon), then it refers to all the entries in that dimension of the matrix.

```{code-cell}
A(:, 1:2:end)        % all of the odd columns
```

:::{index} ! MATLAB; diag
:::

The matrix and vector senses of addition, subtraction, scalar multiplication, multiplication, and power are all handled by the usual symbols. 
```{tip}
:class: dropdown
Use `diag` to construct a matrix by its diagonals. A more general syntax puts elements on super- or subdiagonals.
```

```{code-cell}
B = diag([-1, 0, -5])   % create a diagonal matrix
```

```{code-cell}
size(A)
size(B)
```

```{code-cell}
BA = B * A     % matrix product
```

`A * B` causes an error here, because the dimensions aren't compatible.
```{tip}
:class: dropdown
Errors are formally called *exceptions* in Julia.
```

```{code-cell} julia
:tags: raises-exception
A * B    % throws an error
```

A square matrix raised to an integer power is the same as repeated matrix multiplication.

```{code-cell}
B^3    % same as B*B*B
```

Sometimes one instead wants to treat a matrix or vector as a mere array and simply apply a single operation to each element of it. For multiplication, division, and power, the corresponding operators start with a dot.

```{code-cell}
C = -A;
```

Because both matrices are $3\times 5$, `A * C` would be an error here, but elementwise operations are fine.

```{code-cell}
elementwise = A .* C
```

```{index} MATLAB; broadcasting
```

The two operands of a dot operator have to have the same size—unless one is a scalar, in which case it is expanded or *broadcast* to be the same size as the other operand.

```{code-cell}
x_to_two = x .^ 2
```

```{code-cell}
two_to_x = 2 .^ x
```

```{tip}
:class: dropdown
Most of the mathematical functions, such as cos, sin, log, exp, and sqrt, can operate elementwise on vectors and matrices. 
```

```{code-cell}
cos(pi * x) 
```
``````

### 2.3 @section-linsys-linear-systems
(demo-systems-backslash-matlab)=
``````{dropdown} @demo-systems-backslash
For a square matrix $\mathbf{A}$, the syntax `A \ b` is mathematically equivalent to $\mathbf{A}^{-1} \mathbf{b}$. 

```{code-cell}
A = [1 0 -1; 2 2 1; -1 -3 0]
```

```{code-cell}
b = [1; 2; 3]
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

If the matrix $\mathbf{A}$ is singular, you may get a warning and nonsense result.

```{code-cell} julia
:tags: raises-exception
A = [0 1; 0 0]
b = [1; -1]
x = A \ b
```

In this case, we can check that the rank of $\mathbf{A}$ is less than its number of columns, indicating singularity.
```{tip}
:class: dropdown
The function `rank` computes the rank of a matrix. However, it is numerically unstable for matrices that are nearly singular, in a sense to be defined in a later section.
```

```{code-cell}
rank(A)
```

A linear system with a singular matrix might have no solution or infinitely many solutions, but in either case, backslash will fail. Moreover, detecting singularity is a lot like checking whether two floating-point numbers are *exactly* equal: because of roundoff, it could be missed. In {numref}`section-linsys-condition-number` we'll find a robust way to fully describe this situation.
``````

(demo-systems-triangular-matlab)=
``````{dropdown} @demo-systems-triangular
```{index} ! MATLAB; tril, ! MATLAB; triu
```

It's easy to get just the lower triangular part of any matrix using the `tril` function.
```{tip}
:class: dropdown
Use `tril` to return a matrix that zeros out everything above the main diagonal. The `triu` function zeros out below the diagonal.
```

```{code-cell}
A = randi(9, 5, 5);
L = tril(A)
```

We'll set up and solve a linear system with this matrix.

```{code-cell}
b = ones(5);
x = forwardsub(L, b)
```

```{index} residual
```

It's not clear how accurate this answer is. However, the residual should be zero or comparable to $\macheps$.

```{code-cell}
b - L * x
```

```{index} ! MATLAB; diag, ! MATLAB; eye
```

Next, we'll engineer a problem to which we know the exact answer. 
```{tip}
:class: dropdown
The `eye` function creates an identity matrix. The `diag` function uses 0 as the main diagonal, positive integers as superdiagonals, and negative integers as subdiagonals.
```

```{code-cell}
alpha = 0.3;
beta = 2.2;
U = eye(5) + diag([-1 -1 -1 -1], 1);
U(1, [4, 5]) = [alpha - beta, beta]
```

```{code-cell}
x_exact = ones(5);
b = [alpha; 0; 0; 0; 1];
```

Now we use backward substitution to solve for $\mathbf{x}$, and compare to the exact solution we know already.

```{code-cell}
x = backsub(U, b);
err = x - x_exact
```

Everything seems OK here. But another example, with a different value for $\beta$, is more troubling.

```{code-cell}
alpha = 0.3;
beta = 1e12;
U = eye(5) + diag([-1 -1 -1 -1], 1);
U(1, [4, 5]) = [alpha - beta, beta];
b = [alpha; 0; 0; 0; 1];

x = backsub(U, b);
err = x - x_exact
```

It's not so good to get 4 digits of accuracy after starting with sixteen! The source of the error is not hard to track down. Solving for $x_1$ performs $(\alpha-\beta)+\beta$ in the first row. Since $|\alpha|$ is so much smaller than $|\beta|$, this a recipe for losing digits to subtractive cancellation.
``````
### 2.4 @section-linsys-lu
(demo-lu-outertri-matlab)= 
``````{dropdown} @demo-lu-outertri
:open: false

```{index} MATLAB; tril, MATLAB; triu
```
We explore the outer product formula for two random triangular matrices.

```{code-cell}
L = tril( randi(9, 3, 3) )
```

```{code-cell}
U = triu( randi(9, 3, 3) )
```

Here are the three outer products in the sum in {eq}`matrixouter`:

```{code-cell}
L(:, 1) * U(1, :)
```

```{code-cell}
L(:, 2) * U(2, :)
```

```{code-cell}
L(:, 3) * U(3, :)
```

Simply because of the triangular zero structures, only the first outer product contributes to the first row and first column of the entire product. 
``````

(demo-lu-derive-matlab)=
``````{dropdown} @demo-lu-derive
For illustration, we work on a $4 \times 4$ matrix. We name it with a subscript in preparation for what comes.

```{code-cell}
A_1 = [
     2    0    4     3 
    -4    5   -7   -10 
     1   15    2   -4.5
    -2    0    2   -13
    ];
L = eye(4);
U = zeros(4, 4);
```

Now we appeal to {eq}`outer-row1`. Since $L_{11}=1$, we see that the first row of $\mathbf{U}$ is just the first row of $\mathbf{A}_1$.

```{code-cell}
U(1, :) = A_1(1, :)
```

From {eq}`outer-col1`, we see that we can find the first column of $\mathbf{L}$ from the first column of $\mathbf{A}_1$. 

```{code-cell}
L(:, 1) = A_1(:, 1) / U(1, 1)
```

 We have obtained the first term in the sum {eq}`matrixouter` for $\mathbf{L}\mathbf{U}$, and we subtract it away from $\mathbf{A}_1$.

```{code-cell}
A_2 = A_1 - L(:, 1) * U(1, :)
```

Now $\mathbf{A}_2 = \boldsymbol{\ell}_2\mathbf{u}_2^T + \boldsymbol{\ell}_3\mathbf{u}_3^T + \boldsymbol{\ell}_4\mathbf{u}_4^T.$ If we ignore the first row and first column of the matrices in this equation, then in what remains we are in the same situation as at the start. Specifically, only $\boldsymbol{\ell}_2\mathbf{u}_2^T$ has any effect on the second row and column, so we can deduce them now.

```{code-cell}
U(2, :) = A_2(2, :)
L(:, 2) = A_2(:, 2) / U(2, 2)
```

If we subtract off the latest outer product, we have a matrix that is zero in the first *two* rows and columns. 

```{code-cell}
A_3 = A_2 - L(:, 2) * U(2, :)
```

Now we can deal with the lower right $2\times 2$ submatrix of the remainder in a similar fashion.

```{code-cell}
U(3, :) = A_3(3, :);
L(:, 3) = A_3(:, 3) / U(3, 3);
A_4 = A_3 - L(:, 3) * U(3, :)
```

Finally, we pick up the last unknown in the factors.

```{code-cell}
U(4, 4) = A_4(4, 4);
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
A_1 - L * U
```

In floating point, we cannot expect the difference to be exactly zero as we found in this toy example. Instead, we would be satisfied to see that each element of the difference above is comparable in size to machine precision.

``````

(demo-lu-solve-matlab)=
``````{dropdown} @demo-lu-solve
Here are the data for a linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$. 

```{code-cell}
A = [2 0 4 3; -4 5 -7 -10; 1 15 2 -4.5; -2 0 2 -13];
b = [4; 9; 9; 4];
```

We apply {numref}`Function {number} <function-lufact>` and then do two triangular solves.

```{code-cell}
[L, U] = lufact(A)
z = forwardsub(L, b);
x = backsub(U, z);
```

A check on the residual assures us that we found the solution.

```{code-cell}
b - A * x
```
``````

### 2.5 @section-linsys-efficiency
(demo-flops-mvmult-matlab)=
``````{dropdown} @demo-flops-mvmult
Here is a straightforward implementation of matrix-vector multiplication.

```{code-cell}
n = 6;
A = magic(n);
x = ones(n,1);
y = zeros(n,1);
for i = 1:n
    for j = 1:n
        y(i) = y(i) + A(i,j)*x(j);   % 2 flops
    end
end
```

Each of the loops implies a summation of flops. The total flop count for this algorithm is

$$ \sum_{i=1}^n \sum_{j=1}^n 2 = \sum_{i=1}^n 2n = 2n^2. $$

Since the matrix $\mathbf{A}$ has $n^2$ elements, all of which have to be involved in the product, it seems unlikely that we could get a flop count that is smaller than $O(n^2)$ in general.

```{index} ! MATLAB; tic and toc
```

Let's run an experiment with the built-in matrix-vector multiplication, using `tic` and `toc` to time the operation.

```{code-cell}
n_ = (400:400:4000)';
t_ = zeros(size(n_));
for i = 1:length(n_)
    n = n_(i);
    A = randn(n, n);  x = randn(n, 1);
    tic    % start a timer
    for j = 1:100      % repeat 100 times
        A*x;
    end
    t = toc;           % read the timer
    t_(i) = t / 100;   % seconds per instance
end
```

The reason for doing multiple repetitions at each value of $n$ in the loop above is to avoid having times so short that the resolution of the timer is significant.

```{code-cell}
table(n_, t_, 'variablenames', {'size', 'time'})
```

```{index} MATLAB; Boolean indexing
```

Looking at the timings just for $n=2000$ and $n=4000$, they have ratio
```{tip}
:class: dropdown
The expression `n_==4000` here produces a vector of Boolean (true/false) values the same size as `n_`. This result is used to index within `t_`, accessing only the value for which the comparison is true.
```

```{code-cell}
t_(n_==4000) / t_(n_==2000)
```

If the run time is dominated by flops, then we expect this ratio to be 

$$
\frac{2(4000)^2}{2(2000)^2}=4.
$$
``````

(demo-flops-loglog-matlab)=
``````{dropdown} @demo-flops-loglog

Let's repeat the previous experiment for more, and larger, values of $n$.

```{code-cell}
n_ = (400:400:6000)';
t_ = zeros(size(n_));
for i = 1:length(n_)
    n = n_(i);
    A = randn(n, n);  x = randn(n, 1);
    tic    % start a timer
    for j = 1:100      % repeat ten times
        A*x;
    end
    t = toc;          % read the timer
    t_(i) = t / 100;   % seconds per instance
end
```

Plotting the time as a function of $n$ on log-log scales is equivalent to plotting the logs of the variables.

```{code-cell}
clf    % clear any existing figure
loglog(n_, t_, '.-')
xlabel('size of matrix')
ylabel('time (sec)')
title(('Timing of matrix-vector multiplications'));
```

You can see that while the full story is complicated, the graph is trending to a straight line of positive slope. For comparison, we can plot a line that represents $O(n^2)$ growth exactly. (All such lines have slope equal to 2.)

```{code-cell}
hold on
loglog(n_, t_(1) * (n_ / n_(1)).^2, '--')
axis tight
legend('data', 'O(n^2)', 'location', 'southeast');
```
``````

(demo-flops-lufact-matlab)=
``````{dropdown} @demo-flops-lufact
We'll test the conclusion of $O(n^3)$ flops experimentally, using the built-in `lu` function instead of the purely instructive `lufact`.
```{tip}
:class: dropdown
The first time a function is invoked, there may be significant time needed to compile it in memory. Thus, when timing a function, run it at least once before beginning the timing.
```

```{code-cell}
n_ = (200:100:2400)';
t_ = zeros(size(n_));
for i = 1:length(n_)
    n = n_(i);
    A = randn(n, n);  
    tic    % start a timer
    for j = 1:6,  [L, U] = lu(A);  end
    t = toc;
    t_(i) = t / 6;  
end
```

We plot the timings on a log-log graph and compare it to $O(n^3)$. The result could vary significantly from machine to machine, but in theory the data should start to parallel the line as $n\to\infty$.

```{code-cell}
clf
loglog(n_,t_,'.-')
hold on, loglog(n_,t_(end)*(n_/n_(end)).^3,'--')
axis tight
xlabel('size of matrix'), ylabel('time (sec)')
title('Timing of LU factorization')
legend('lu','O(n^3)','location','southeast');
```
``````

### 2.6 @section-linsys-pivoting
(demo-pivoting-fail-matlab)=
``````{dropdown} @demo-pivoting-fail
Here is a previously encountered matrix that factors well.

```{code-cell}
A = [
    2 0 4 3
    -4 5 -7 -10
    1 15 2 -4.5
    -2 0 2 -13
    ];
[L, U] = lufact(A);
L
```

If we swap the second and fourth rows of $\mathbf{A}$, the result is still nonsingular. However, the factorization now fails.

```{code-cell}
A([2, 4], :) = A([4, 2], :);    % swap rows 2 and 4
[L, U] = lufact(A);
L
```

```{index} MATLAB; NaN
```

The presence of `NaN` in the result indicates that some impossible operation was required. The source of the problem is easy to locate. We can find the first outer product in the factorization just fine:

```{code-cell}
U(1, :) = A(1, :);
L(:, 1) = A(:, 1) / U(1, 1)
A = A - L(:, 1) * U(1, :)
```

The next step is `U(2, :) = A(2, :)`, which is also OK. But then we are supposed to divide by `U(2, 2)`, which is zero. The algorithm cannot continue.
``````

(demo-pivoting-fix-matlab)=
``````{dropdown} @demo-pivoting-fix
Here is the trouble-making matrix from {numref}`Demo {number} <demo-pivoting-fail>`.

```{code-cell}
A_1 = [2 0 4 3; -2 0 2 -13; 1 15 2 -4.5; -4 5 -7 -10]
```

```{index} ! MATLAB; max, ! MATLAB; \~
```

We now find the largest candidate pivot in the first column. We don't care about sign, so we take absolute values before finding the max.
```{tip}
:class: dropdown
The second output of `max` returns the location of the largest element of a vector. The `~` symbol is used to ignore the value of the first output.
```


```{code-cell}
[~, i] = max( abs(A_1(:, 1)) ) 
```

This is the row of the matrix that we extract to put into $\mathbf{U}$. That guarantees that the division used to find $\boldsymbol{\ell}_1$ will be valid.

```{code-cell}
L = zeros(4, 4);
U = zeros(4, 4);
U(1, :) = A_1(i, :);
L(:, 1) = A_1(:, 1) / U(1, 1);
A_2 = A_1 - L(:, 1) * U(1, :)
```

Observe that $\mathbf{A}_2$ has a new zero row and zero column, but the zero row is the fourth rather than the first. However, we forge on by using the largest possible pivot in column 2 for the next outer product.

```{code-cell}
[~, i] = max( abs(A_2(:, 2)) )
U(2, :) = A_2(i, :);
L(:, 2) = A_2(:, 2) / U(2, 2);
A_3 = A_2 - L(:, 2) * U(2, :)
```

Now we have zeroed out the third row as well as the second column. We can finish out the procedure.

```{code-cell}
[~, i] = max( abs(A_3(:, 3)) ) 
U(3, :) = A_3(i, :);
L(:, 3) = A_3(:, 3) / U(3, 3);
A_4 = A_3 - L(:, 3) * U(3, :)
```

```{code-cell}
[~, i] = max( abs(A_4(:, 4)) ) 
U(4, :) = A_4(i, :);
L(:, 4) = A_4(:, 4) / U(4, 4);
```

We do have a factorization of the original matrix:

```{code-cell}
A_1 - L * U
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

(demo-pivoting-permute-matlab)=
``````{dropdown} @demo-pivoting-permute
Here again is the matrix from {numref}`Demo {number} <demo-pivoting-fix>`.

```{code-cell}
A = [2 0 4 3; -2 0 2 -13; 1 15 2 -4.5; -4 5 -7 -10]
```

As the factorization proceeded, the pivots were selected from rows 4, 3, 2, and finally 1. If we were to put the rows of $\mathbf{A}$ into that order, then the algorithm would run exactly like the plain LU factorization from {numref}`section-linsys-lu`. 

```{code-cell}
B = A([4, 3, 2, 1], :);
[L, U] = lufact(B);
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

(demo-pivoting-usage-matlab)=
``````{dropdown} @demo-pivoting-usage
The third output of `plufact` is the permutation vector we need to apply to $\mathbf{A}$.

```{code-cell}
A = randi(20, 4, 4);
[L, U, p] = plufact(A);
A(p, :) - L * U    % should be ≈ 0
```

Given a vector $\mathbf{b}$, we solve $\mathbf{A}\mathbf{x}=\mathbf{b}$ by first permuting the entries of $\mathbf{b}$ and then proceeding as before.

```{code-cell}
b = rand(4, 1);
z = forwardsub(L, b(p));
x = backsub(U, z)
```

A residual check is successful:

```{code-cell}
b - A*x
```
``````

(demo-pivoting-builtin-matlab)=
``````{dropdown} @demo-pivoting-builtin
With the syntax `A \ b`, the matrix `A` is PLU-factored, followed by two triangular solves.

```{code-cell}
A = randn(500, 500);    % 500x500 with normal random entries
tic; for k=1:50; A \ rand(500, 1); end; toc
```

In {numref}`section-linsys-efficiency` we showed that the factorization is by far the most costly part of the solution process. A factorization object allows us to do that costly step only once per unique matrix. 

```{code-cell}
[L, U, p] = lu(A, 'vector');    % keep factorization result
tic
for k=1:50
    b = rand(500, 1);
    U \ (L \ b(p));
end
toc
```
``````

(demo-pivoting-stable-matlab)=
``````{dropdown} @demo-pivoting-stable
We construct a linear system for this matrix with $\epsilon=10^{-12}$ and exact solution $[1, 1]$:

```{code-cell}
ep = 1e-12
A = [-ep 1; 1 -1];
b = A * [1; 1];
```

We can factor the matrix without pivoting and solve for $\mathbf{x}$.

```{code-cell}
[L, U] = lufact(A);
x = backsub( U, forwardsub(L, b) )
```

Note that we have obtained only about 5 accurate digits for $x_1$. We could make the result even more inaccurate by making $\epsilon$ even smaller:

```{code-cell}
ep = 1e-20; A = [-ep 1; 1 -1];
b = A * [1; 1];
[L, U] = lufact(A);
x = backsub( U, forwardsub(L, b) )
```

This effect is not due to ill conditioning of the problem—a solution with PLU factorization works perfectly:

```{code-cell}
A \ b
```
``````

### 2.7 @section-linsys-norms
(demo-norms-vector-matlab)=
``````{dropdown} @demo-norms-vector
```{index} ! MATLAB; norm
```

```{code-cell}
x = [2; -3; 1; -1];
twonorm = norm(x)    % or norm(x, 2)
```

```{code-cell}
infnorm = norm(x, Inf)
```

```{code-cell}
onenorm = norm(x, 1)
```
``````

(demo-norms-matrix-matlab)=
``````{dropdown} @demo-norms-matrix

```{code-cell}
A = [ 2 0; 1 -1 ]
```

```{index} ! MATLAB; norm
```

The default matrix norm is the 2-norm.

```{code-cell}
twonorm = norm(A)
```

You can get the 1-norm as well.

```{code-cell}
onenorm = norm(A, 1)
```

```{index} ! MATLAB; max, ! MATLAB; sum
```

According to {eq}`mxonenorm`, the matrix 1-norm is equivalent to the maximum of the sums down the columns (in absolute value).
```{tip}
:class: dropdown
Use `sum` to sum along a dimension of a matrix. The `max` and `min` functions also work along one dimension.
```

```{code-cell}
% Sum down the rows (1st matrix dimension):
max( sum(abs(A), 1) )   
```

Similarly, we can get the $\infty$-norm and check our formula for it.

```{code-cell}
infnorm = norm(A, Inf)
```

```{code-cell}
% Sum across columns (2nd matrix dimension):
max( sum(abs(A), 2) )  
```
Next we illustrate a geometric interpretation of the 2-norm. First, we will sample a lot of vectors on the unit circle in $\mathbb{R}^2$.
```{tip}
:class: dropdown
You can use functions as values, e.g., as elements of a vector. 
```

```{index} ! MATLAB; subplot
```

```{code-cell}
theta = linspace(0, 2*pi, 601);
x = [ cos(theta); sin(theta) ];    % 601 unit column vectors
clf
subplot(1, 2, 1)
plot(x(1, :), x(2, :)), axis equal
title('Unit circle in 2-norm')
xlabel('x_1')
ylabel(('x_2'));
```

The linear function $\mathbf{f}(\mathbf{x}) = \mathbf{A}\mathbf{x}$ defines a mapping from $\mathbb{R}^2$ to $\mathbb{R}^2$. We can apply `A` to every column of `x` by using a single matrix multiplication.

```{code-cell}
Ax = A * x;
```

The image of the transformed vectors is an ellipse that just touches the circle of radius $\|\mathbf{A}\|_2$:

```{code-cell}
subplot(1,2,2), plot(Ax(1,:), Ax(2,:)), axis equal
hold on, plot(twonorm * x(1,:), twonorm * x(2,:), '--')
title('Image of Ax, with ||A||')
xlabel('x_1')
ylabel(('x_2'));
```
``````


### 2.8 @section-linsys-condition-number
(demo-condition-bound-matlab)=
``````{dropdown} @demo-condition-bound

```{index} ! MATLAB; cond
```
MATLAB has a function `cond` to compute matrix condition numbers. By default, the 2-norm is used. As an example, the family of *Hilbert matrices* is famously badly conditioned. Here is the $6\times 6$ case.

```{code-cell}
A = hilb(6)
kappa = cond(A)
```

Because $\kappa\approx 10^8$, it's possible to lose nearly 8 digits of accuracy in the process of passing from $\mathbf{A}$ and $\mathbf{b}$ to $\mathbf{x}$. That fact is independent of the algorithm; it's inevitable once the data are expressed in finite precision. 

Let's engineer a linear system problem to observe the effect of a perturbation. We will make sure we know the exact answer.

```{code-cell}
x = (1:6)';
b = A * x;
```

Now we perturb the system matrix and vector randomly by $10^{-10}$ in norm.

```{code-cell} 
dA = randn(size(A));  dA = 1e-10 * (dA / norm(dA));
db = randn(size(b));  db = 1e-10 * (db / norm(db));
```

We solve the perturbed problem using pivoted LU and see how the solution was changed.

```{code-cell}
new_x = ((A + dA) \ (b + db));
dx = new_x - x;
```

Here is the relative error in the solution.

```{code-cell}
relative_error = norm(dx) / norm(x)
```

And here are upper bounds predicted using the condition number of the original matrix.

```{code-cell}
upper_bound_b = (kappa * norm(db) / norm(b))
upper_bound_A = (kappa * norm(dA) / norm(A))
```

Even if we didn't make any manual perturbations to the data, machine roundoff does so at the relative level of $\macheps$.

```{code-cell}
dx = A\b - x;
relative_error = norm(dx) / norm(x)
rounding_bound = kappa * eps
```

Larger Hilbert matrices are even more poorly conditioned:

```{code-cell}
A = hilb(14);
kappa = cond(A)
```

Note that $\kappa$ exceeds $1/\macheps$. In principle we therefore may end up with an answer that has relative error greater than 100%.

```{code-cell}
rounding_bound = kappa * eps
```

Let's put that prediction to the test.

```{code-cell}
x = (1:14)';  b = A * x;
dx = A\b - x;
relative_error = norm(dx) / norm(x)
```

As anticipated, the solution has zero accurate digits in the 2-norm.
``````

### 2.9 @section-linsys-structure
(demo-structure-banded-matlab)=
``````{dropdown} @demo-structure-banded
```{index} ! MATLAB; fill, MATLAB; diagm, ! MATLAB; diag
```

Here is a small tridiagonal matrix. Note that there is one fewer element on the super- and subdiagonals than on the main diagonal.

```{code-cell}
A = [ 2 -1  0  0  0  0
      4  2 -1  0  0  0
      0  3  0 -1  0  0
      0  0  2  2 -1  0
      0  0  0  1  1 -1
      0  0  0  0  0  2 ];
```

```{index} ! MATLAB; diag
```

We can extract the elements on any diagonal using the `diag` function. The main or central diagonal is numbered zero, above and to the right of that is positive, and below and to the left is negative.
```{tip}
:class: dropdown
The `diag` function extracts the elements from a specified diagonal of a matrix.
```

```{code-cell}
diag_main = diag(A, 0)'
diag_plusone = diag(A, 1)'
diag_minusone = diag(A,-1)'
```

% We can also put whatever numbers we like onto any diagonal with `diag`.
    
```{code-cell}
A = A + diag([5 8 6 7], 2)
```

The lower and upper bandwidths of $\mathbf{A}$ are repeated in the factors from the unpivoted LU factorization. 

```{code-cell}
[L, U] = lufact(A)
```
``````

(demo-structure-symm-matlab)=
``````{dropdown} @demo-structure-symm

We begin with a symmetric $\mathbf{A}$.

```{code-cell}
A_1 = [ 2     4     4     2
        4     5     8    -5
        4     8     6     2
        2    -5     2   -26 ];
```

We won't use pivoting, so the pivot element is at position (1,1). This will become the first element on the diagonal of $\mathbf{D}$. Then we divide by that pivot to get the first column of $\mathbf{L}$.

```{code-cell}
L = eye(4);
d = zeros(4, 1);
d(1) = A_1(1, 1);
L(:, 1) = A_1(:, 1) / d(1);
A_2 = A_1 - d(1) * L(:, 1) * L(:, 1)'
```

We are now set up the same way for the submatrix in rows and columns 2–4.

```{code-cell}
d(2) = A_2(2, 2);
L(:, 2) = A_2(:, 2) / d(2);
A_3 = A_2 - d(2) * L(:, 2) * L(:, 2)'
```

We continue working our way down the diagonal.

```{code-cell}
d(3) = A_3(3, 3);
L(:, 3) = A_3(:, 3) / d(3);
A_4 = A_3 - d(3) * L(:, 3) * L(:, 3)'
d(4) = A_4(4, 4);
d
L
```

We have arrived at the desired factorization, which we can validate:

```{code-cell}
norm(A_1 - (L * diag(d) * L'))
```
``````

(demo-structure-cholesky-matlab)=
``````{dropdown} @demo-structure-cholesky
A randomly chosen matrix is extremely unlikely to be symmetric. However, there is a simple way to symmetrize one.

```{code-cell}
A = magic(4) + eye(4);
B = A + A'
```

```{index} ! MATLAB; chol
```

Similarly, a random symmetric matrix is unlikely to be positive definite. The Cholesky algorithm always detects a non-PD matrix by quitting with an error.
```{tip}
:class: dropdown
The `chol` function computes a Cholesky factorization if possible, or throws an error for a non-positive-definite matrix. 
```

```{warning} 
The `chol` function does *not* check for symmetry. It may give a nonsensical result if the input is not symmetric.
```

```{code-cell}
:tags: raises-exception
chol(B)    % throws an error
```

It's not hard to manufacture an SPD matrix to try out the Cholesky factorization.

```{code-cell}
B = A' * A;
R = chol(B)
```

Here we validate the factorization:

```{code-cell}
norm(R' * R - B) / norm(B)
```
``````