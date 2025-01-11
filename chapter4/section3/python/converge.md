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
[**Demo %s**](#demo-newton-converge)

We again look at finding a solution of $x e^x=2$ near $x=1$. To apply Newton's method, we need to calculate values of both the residual function $f$ and its derivative.

```{code-cell}
f = lambda x: x * exp(x) - 2
df_dx = lambda x: exp(x) * (x + 1)
```

We don't know the exact root, so we use `nlsolve` to determine a proxy for it.

```{code-cell}
r = root_scalar(f, bracket=[0.8, 1.0]).root
print(r)
```

We use $x_1=1$ as a starting guess and apply the iteration in a loop, storing the sequence of iterates in a vector.

```{code-cell}
x = ones(5)
for k in range(4):
    x[k + 1] = x[k] - f(x[k]) / df_dx(x[k])

print(x)
```

Here is the sequence of errors.

```{code-cell}
err = x - r
print(err)
```

The exponents in the scientific notation definitely suggest a squaring sequence. We can check the evolution of the ratio in {eq}`quadratictest`.

```{code-cell}
logerr = log(abs(err))
for i in range(len(err) - 1):
    print(logerr[i+1] / logerr[i])
```

The clear convergence to 2 above constitutes good evidence of quadratic convergence.
