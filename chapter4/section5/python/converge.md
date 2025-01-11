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
[**Demo %s**](#demo-newtonsys-converge)

A system of nonlinear equations is defined by its residual and Jacobian.

```{code-cell}
def func(x):
    return array([
        exp(x[1] - x[0]) - 2, 
        x[0] * x[1] + x[2], 
        x[1] * x[2] + x[0]**2 - x[1]
    ])

def jac(x):
    return array([
            [-exp(x[1] - x[0]), exp(x[1] - x[0]), 0],
            [x[1], x[0], 1],
            [2 * x[0], x[2] - 1, x[1]],
    ])
```

Our initial guess at a root is the origin. 

```{code-cell}
x1 = zeros(3)
x = FNC.newtonsys(func, jac, x1)
print(x)
```

The output has one row per iteration, so the last row contains the final Newton estimate. Let's compute its residual.

```{code-cell}
r = x[-1]
f = func(r)
print("final residual:", f)
```

Let's check the convergence rate:

```{code-cell}
logerr = [log(norm(x[k] - r)) for k in range(x.shape[0] - 1)]
for k in range(len(logerr) - 1):
    print(logerr[k+1] / logerr[k])
```

The ratio is apparently converging toward 2, as expected for quadratic convergence.
