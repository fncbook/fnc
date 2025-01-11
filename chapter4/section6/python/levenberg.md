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
[**Demo %s**](#demo-quasi-levenberg)

To solve a nonlinear system, we need to code only the function defining the system, and not its Jacobian.

```{code-cell}
def func(x):
    return array([
        exp(x[1] - x[0]) - 2, 
        x[0] * x[1] + x[2], 
        x[1] * x[2] + x[0]**2 - x[1]
    ])
```

In all other respects usage is the same as for the `newtonsys` function.

```{code-cell}
x1 = zeros(3)
x = FNC.levenberg(func, x1)
print(f"Took {len(x) - 1} iterations.")
```

It's always a good idea to check the accuracy of the root, by measuring the residual (backward error).

```{code-cell}
r = x[-1]
print("backward error:", norm(func(r)))
```
Looking at the convergence in norm, we find a convergence rate between linear and quadratic, like with the secant method:

```{code-cell}
logerr = [log(norm(x[k] - r)) for k in range(len(x) - 1)]
for k in range(len(logerr) - 1):
    print(logerr[k+1] / logerr[k])
```
