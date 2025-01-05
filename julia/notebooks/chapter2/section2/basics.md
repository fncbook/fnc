---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
:::{index} ! Julia; size, ! Julia; length
:::

In Julia, vectors and matrices are one-dimensional and two-dimensional arrays, respectively. Square brackets are used to enclose elements of a matrix or vector. Use spaces for horizontal concatenation, and semicolons or new lines to indicate vertical concatenation.
```{tip}
The `size` function returns the number of rows and columns in a matrix. Use `length` to get the number of elements in a vector or matrix.
```

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

Accessing an element is done by giving one (for a vector) or two (for a matrix) index values within square brackets. 
```{tip}
The `end` keyword refers to the last element in a dimension. It saves you from having to compute and store the size of the matrix first.
```
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
`A * B` causes an error here, because the dimensions aren't compatible.
```{tip}
Errors are formally called *exceptions* in Julia.
```

```{code-cell} julia
:tags: [raises-exception]
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
x_to_two = x .^ 2
```

```{code-cell}
two_to_x = 2 .^ x
```

Most of the mathematical functions, such as cos, sin, log, exp, and sqrt, expect scalars as operands. However, you can broadcast any function, including ones that you have defined, across a vector or array by using a special dot syntax.
```{tip}
A dot added to the end of a function name means to apply the function elementwise to an array.
```

```{code-cell}
show(cos.(π * x))    # broadcast to a function
```

Rather than dotting multiple individual functions, you can use `@.` before an expression to broadcast everything within it.

```{code-cell}
show(@. cospi( (x + 1)^3) )    # broadcast an entire expression
```
