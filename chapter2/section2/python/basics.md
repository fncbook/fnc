---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
exec(open("../../../python/FNC_init.py").read())
```
[**Demo %s**](#demo-matrices-basics)

```{note}
While NumPy does have distinct representations for matrices and 2D arrays, use of the explicit matrix class is officially discouraged. We follow this advice here and use arrays to represent both matrices and vectors.
```

:::{index} ! Python; array, ! Python; shape
:::

<!-- :::{index}
see: Python; size, Python; shape
::: -->

A vector is created using square brackets and commas to enclose and separate its entries.

```{code-cell} 
x = array([3, 3, 0, 1, 0 ])
print(x.shape)
```

To construct a matrix, you nest the brackets to create a "vector of vectors". The inner vectors are the rows.

```{code-cell} 
A = array([ 
    [1, 2, 3, 4, 5],
    [50, 40, 30, 20, 10], 
    [pi, sqrt(2), exp(1), (1+sqrt(5))/2, log(3)] 
    ])

print(A)
print(A.shape)
```

In this text, we treat all vectors as equivalent to matrices with a single column. That isn't true in NumPy, because even an $n \times 1$ array has two dimensions, unlike a vector.

```{code-cell} 
array([[3], [1], [2]]).shape
```

:::{index} ! Python; hstack, ! Python; vstack
:::

You can concatenate arrays with compatible dimensions using `hstack` and `vstack`.

```{code-cell} 
print( hstack([A, A]) )
```

```{code-cell} 
print( vstack([A, A]) )
```

```{index} ! Python; transpose, ! Python; adjoint
```

Transposing a matrix is done by appending `.T` to it. 

```{code-cell} 
print(A.T)
```

For matrices with complex values, we usually want instead the adjoint or hermitian, which is `.conj().T`. 

```{code-cell} 
print((x + 1j).conj().T)
```

:::{index} ! Python; arange, ! Python; linspace
:::

There are many convenient shorthand ways of building vectors and matrices other than entering all of their entries directly or in a loop. To get a vector with evenly spaced entries between two endpoints, you have two options.

```{code-cell} 
print(arange(1, 7, 2))   # from 1 to 7 (not inclusive), step by 2        
```

```{code-cell} 
print(linspace(-1, 1, 5))   # from -1 to 1 (inclusive), with 5 total values
```

The practical difference between these is whether you want to specify the step size in `arange` or the number of points in `linspace`.

Accessing an element is done by giving one (for a vector) or two index values in square brackets. **In Python, indexing always starts with zero, not 1.**

```{code-cell} 
A = array([ 
    [1, 2, 3, 4, 5],
    [50, 40, 30, 20, 10], 
    linspace(-5, 5, 5) 
    ])
x = array([3, 2, 0, 1, -1 ])
```

```{code-cell} 
print("row 2, col 3 of A:", A[1, 2])
print("first element of x:", x[0])
```

```{index} ! Python; slice, ! Python; \:
```
:::{index} ! Python; indexing arrays
:::

The indices can be ranges, in which case a **slice** or block of the matrix is accessed. You build these using a colon in the form `start:stop`. However, the last value of this range is `stop-1`, not `stop`.

```{code-cell} 
print(A[1:3, 0:2])    # rows 2 and 3, cols 1 and 2
```

If `start` or `stop` is omitted, the range extends to the first or last index.

```{code-cell} 
print(x[1:])  # elements 2 through the end
```

```{code-cell} 
print(A[:2, 0])  # first two rows in column 1
```

Notice in the last case above that even when the slice is in the shape of a column vector, the result is just a vector with one dimension and neither row nor column shape.

There are more variations on the colon ranges. A negative value means to count from the end rather than the beginning. And a colon by itself means to include everything from the relevant dimension.

```{code-cell} 
print(A[:-1, :])    # all rows up to the last, all columns
```

Finally, `start:stop:step` means to step size or stride other than one. You can mix this with the other variations.

```{code-cell} 
print(x[::2])  # all the odd indexes
```

```{code-cell} 
print(A[:, ::-1])  # reverse the columns
```

The matrix and vector senses of addition, subtraction, and scalar multiplication and division are all handled by the usual symbols. Two matrices of the same size (what NumPy calls shape) are operated on elementwise. 

```{code-cell} 
print(A - 2 * ones([3, 5]))  # subtract two from each element
```

```{index} ! Python; broadcast
```

If one operand has a smaller number of dimensions than the other, Python tries to **broadcast** it in the "missing" dimension(s), and the operation proceeds if the resulting shapes are identical. 

```{code-cell} 
print(A - 2)    # subtract two from each element
```

```{code-cell} 
u = array([1, 2, 3, 4, 5])
print(A - u)    # repeat this row for every row of A
```

```{code-cell} 
:tags: [raises-exception]
v = array([1, 2, 3])
print(A - v)  # broadcasting this would be 3x3, so it's an error
```

```{code-cell} 
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

```{code-cell} 
B = diag([-1, 0, -5])    # create a diagonal 3x3
print(B @ A)    # matrix product
```

$AB$ is undefined for these matrix sizes. 

```{code-cell} 
:tags: [raises-exception]
print(A @ B)    # incompatible sizes
```

```{index} ! Python; elementwise multiplication, Python; broadcasting
```

The multiplication operator `*` is reserved for elementwise multiplication. Both operands have to be the same size, after any potential broadcasts.

```{code-cell} 
:tags: [raises-exception]
print(B * A)    # not the same size, so it's an error
```

```{code-cell} 
print((A / 2) * A)    # elementwise
```
To raise to a power elementwise, use a double star. This will broadcast as well.

```{code-cell} 
print(B)
print(B**3)
```

```{code-cell} 
print(x)
print(2.0**x)
```

```{danger}
If `A` is a matrix, `A**2` is *not* the same as mathematically raising it to the power 2.
```


```{index} Python; broadcasting
```

Most of the mathematical functions, such as cos, sin, log, exp and sqrt, expecting scalars as operands will be broadcast to arrays.

```{code-cell} 
print(cos(pi * x))      
```
