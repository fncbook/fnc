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
[**Demo %s**](#demo-systems-triangular)


```{index} ! Python; tril, ! Python; triu
```

It's easy to get just the lower triangular part of any matrix using the `tril` function.

```{code-cell} 
A = 1 + floor(9 * random.rand(5, 5))
L = tril(A)
print(L)
```

We'll set up and solve a linear system with this matrix.

```{code-cell} 
b = ones(5)
x = FNC.forwardsub(L, b)
print(x)
```

It's not clear how accurate this answer is. However, the residual should be zero or comparable to $\macheps$.

```{code-cell} 
b - L @ x
```

Next we'll engineer a problem to which we know the exact answer. 

```{code-cell} 
alpha = 0.3;
beta = 2.2;
U = diag(ones(5)) + diag([-1, -1, -1, -1], k=1)
U[0, 3:5] = [ alpha - beta, beta ]
print(U)
```

```{code-cell} 
x_exact = ones(5)
b = array([alpha, 0, 0, 0, 1])
x = FNC.backsub(U, b)
print("error:", x - x_exact)
```

Everything seems OK here. But another example, with a different value for $\beta$, is more troubling.

```{code-cell} 
alpha = 0.3;
beta = 1e12;
U = diag(ones(5)) + diag([-1, -1, -1, -1], k=1)
U[0, 3:5] = [ alpha - beta, beta ]
b = array([alpha, 0, 0, 0, 1])

x = FNC.backsub(U, b)
print("error:", x - x_exact)
```

It's not so good to get 4 digits of accuracy after starting with sixteen! But the source of the error is not hard to track down. Solving for $x_1$ performs $(\alpha-\beta)+\beta$ in the first row. Since $|\alpha|$ is so much smaller than $|\beta|$, this a recipe for losing digits to subtractive cancellation.
