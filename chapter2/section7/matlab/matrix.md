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
[**Demo %s**](#demo-norms-matrix)


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
