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
[**Demo %s**](#demo-condition-bound)


```{index} ! Python; cond
```

The function `cond` from `numpy.linalg` is used to computes matrix condition numbers. By default, the 2-norm is used. As an example, the family of *Hilbert matrices* is famously badly conditioned. Here is the $6\times 6$  case. 

```{code-cell} 
A = array([ 
    [1/(i + j + 2) for j in range(6)] 
    for i in range(6) 
    ])
print(A)
```

```{code-cell} 
from numpy.linalg import cond
kappa = cond(A)
print(f"kappa is {kappa:.3e}")
```

Next we engineer a linear system problem to which we know the exact answer.

```{code-cell} 
x_exact = 1.0 + arange(6)
b = A @ x_exact
```

Now we perturb the data randomly with a vector of norm $10^{-12}$. 

```{code-cell} 
dA = random.randn(6, 6)
dA = 1e-12 * (dA / norm(dA, 2))
db = random.randn(6)
db = 1e-12 * (db / norm(db, 2))
```

We solve the perturbed problem using built-in pivoted LU and see how the solution was changed.

```{code-cell} 
x = linalg.solve(A + dA, b + db) 
dx = x - x_exact
```

Here is the relative error in the solution.

```{code-cell} 
print(f"relative error is {norm(dx) / norm(x_exact):.2e}")
```

And here are upper bounds predicted using the condition number of the original matrix. 

```{code-cell} 
print(f"b_bound: {kappa * 1e-12 / norm(b):.2e}")
print(f"A_bound: {kappa * 1e-12 / norm(A, 2):.2e}")
```

Even if we don't make any manual perturbations to the data, machine epsilon does when we solve the linear system numerically.

```{code-cell} 
x = linalg.solve(A, b)
print(f"relative error: {norm(x - x_exact) / norm(x_exact):.2e}")
print(f"rounding bound: {kappa / 2**52:.2e}")

```

Because $\kappa\approx 10^8$, it's possible to lose 8 digits of accuracy in the process of passing from $A$ and $b$ to $x$. That's independent of the algorithm; it's inevitable once the data are expressed in double precision. 

Larger Hilbert matrices are even more poorly conditioned.

```{code-cell} 
A = array([ [1/(i+j+2) for j in range(14)] for i in range(14) ])
kappa = cond(A)
print(f"kappa is {kappa:.3e}")
```

Before we compute the solution, note that $\kappa$ exceeds `1/eps`. In principle we therefore might end up with an answer that is completely wrong (i.e., a relative error greater than 100%).

```{code-cell} 
print(f"rounding bound: {kappa / 2**52:.2e}")
```

```{code-cell} 
x_exact = 1.0 + arange(14)
b = A @ x_exact  
x = linalg.solve(A, b)
```

We got an answer. But in fact, the error does exceed 100%:

```{code-cell} 
print(f"relative error: {norm(x - x_exact) / norm(x_exact):.2e}")
```
