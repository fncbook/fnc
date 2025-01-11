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
[**Demo %s**](#demo-fp-converge)

We revisit {numref}`Demo %s <demo-fp-spiral>` and investigate the observed convergence more closely. Recall that above we calculated $g'(p)\approx-0.42$ at the convergent fixed point.

```{code-cell}
f = poly1d([1, -4, 3.5])
r = f.roots
print(r)
```

Here is the fixed-point iteration. This time we keep track of the whole sequence of approximations.

```{code-cell}
g = lambda x: x - f(x)
x = zeros(12)
x[0] = 2.1
for k in range(11):
    x[k + 1] = g(x[k])

print(x)
```

It's illuminating to construct and plot the sequence of errors.

```{code-cell}
err = abs(x - max(r))
semilogy(err, "-o")
xlabel("iteration number"), ylabel("error")
title("Convergence of fixed-point iteration");
```

It's quite clear that the convergence quickly settles into a linear rate. We could estimate this rate by doing a least-squares fit to a straight line. Keep in mind that the values for small $k$ should be left out of the computation, as they don't represent the linear trend.

```{code-cell}
p = polyfit(arange(5, 13), log(err[4:]), 1)
print(p)
```

We can exponentiate the slope to get the convergence constant $\sigma$.

```{code-cell}
print("sigma:", exp(p[0]))
```

The error should therefore decrease by a factor of $\sigma$ at each iteration. We can check this easily from the observed data.

```{code-cell}
err[8:] / err[7:-1]
```

The methods for finding $\sigma$ agree well.
