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
[**Demo %s**](#demo-finitediffs-fd-weights)

We will estimate the derivative of $\cos(x^2)$ at $x=0.5$ using five nodes.

```{code-cell}
t = array([0.35, 0.5, 0.57, 0.6, 0.75])   # nodes
f = lambda x: cos(x**2)
dfdx = lambda x: -2 * x * sin(x**2)
exact_value = dfdx(0.5)
```

We have to shift the nodes so that the point of estimation for the derivative is at $x=0$. (To subtract a scalar from a vector, we must use the `.-` operator.)

```{code-cell}
w = FNC.fdweights(t - 0.5, 1)
```

The finite-difference formula is a dot product (i.e., inner product) between the vector of weights and the vector of function values at the nodes.

```{code-cell}
fd_value = dot(w, f(t))
```

We can reproduce the weights in the finite-difference tables by using equally spaced nodes with $h=1$. For example, here is a one-sided formula at four nodes.

```{code-cell}
print(FNC.fdweights(linspace(0, 3, 4), 1))
```
