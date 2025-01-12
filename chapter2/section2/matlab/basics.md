---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
cd  /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init
```
[**Demo %s**](#demo-matrices-basics)

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
:tags: [raises-exception]
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
