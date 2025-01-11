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
[**Demo %s**](#demo-norms-matrix)

```{code-cell} 
from numpy.linalg import norm
A = array([ [2, 0], [1, -1] ])
```

```{index} ! Python; norm
```

The default matrix norm is *not* the 2-norm. Instead, you must provide the 2 explicitly. 

```{code-cell} 
print(norm(A, 2))
```

You can get the 1-norm as well.

```{code-cell} 
print(norm(A, 1))
```

```{index} ! Python; max, Python; sum
```

The 1-norm is equivalent to 

```{code-cell} 
print(max( sum(abs(A), axis=0)) )  # sum down the rows
```

Similarly, we can get the $\infty$-norm and check our formula for it.

```{code-cell} 
print(norm(A, inf))
```

```{code-cell} 
print(max( sum(abs(A), axis=1)) )  # sum across columns 
```

Here we illustrate the geometric interpretation of the 2-norm. First, we will sample a lot of vectors on the unit circle in $\mathbb{R}^2$. 

```{code-cell} 
theta = linspace(0, 2*pi, 601)
x = vstack([cos(theta), sin(theta)])  # 601 unit columns
```

The linear function $\mathbf{f}(\mathbf{x}) = \mathbf{A}\mathbf{x}$ defines a mapping from $\mathbb{R}^2$ to $\mathbb{R}^2$. We can apply `A` to every column of `x` simply by using a matrix multiplication.

```{code-cell} 
y = A @ x
```

We plot the unit circle on the left and the image of all mapped vectors on the right: 

```{code-cell} 
subplot(1,2,1)
plot(x[0, :], x[1, :])
axis("equal")
title("Unit circle")
xlabel("$x_1$")
ylabel("$x_2$")

subplot(1,2,2)
plot(y[0, :], y[1, :])
plot(norm(A, 2) * x[0, :], norm(A,2) * x[1, :],"--")
axis("equal")
title("Image under map")
xlabel("$y_1$")
ylabel("$y_2$");
```

As seen on the right-side plot, the image of the transformed vectors is an ellipse that just touches the circle of radius $\|\mathbf{A}\|_2$.

