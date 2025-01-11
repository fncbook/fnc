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
[**Demo %s**](#demo-fdconverge-order12)

Let's observe the convergence of the formulas in {numref}`Example {number} <example-fd-converge-FD11>` and {numref}`Example {number} <example-fd-converge-FD12>`, applied to the function $\sin(e^{x+1})$ at $x=0$.

```{code-cell}
f = lambda x: sin(exp(x + 1))
exact_value = exp(1) * cos(exp(1))
```

We'll compute the formulas in parallel for a sequence of $h$ values.

```{code-cell}
h_ = array([5 / 10**(n+1) for n in range(6)])
FD = zeros((len(h_), 2))
for (i, h) in enumerate(h_):
    FD[i, 0] = (f(h) - f(0)) / h 
    FD[i, 1] = (f(h) - f(-h)) / (2*h)
results = PrettyTable()
results.add_column("h", h_)
results.add_column("FD1", FD[:, 0])
results.add_column("FD2", FD[:, 1])
print(results)
```

All that's easy to see from this table is that FD2 appears to converge to the same result as FD1, but more rapidly. A table of errors is more informative.

```{code-cell}
errors = FD - exact_value
results = PrettyTable()
results.add_column("h", h_)
results.add_column("error in FD1", errors[:, 0])
results.add_column("error in FD2", errors[:, 1])
print(results)
```

In each row, $h$ is decreased by a factor of 10, so that the error is reduced by a factor of 10 in the first-order method and 100 in the second-order method.

A graphical comparison can be useful as well. On a log-log scale, the error should (as $h\to 0$) be a straight line whose slope is the order of accuracy. However, it's conventional in convergence plots to show $h$ _decreasing_ from left to right, which negates the slopes.

```{code-cell}
plot(h_, abs(errors), "o-", label=["FD1", "FD2"])
gca().invert_xaxis()
# Add lines for perfect 1st and 2nd order.
loglog(h_, h_, "--", label="$O(h)$")
loglog(h_, h_**2, "--", label="$O(h^2)$")
xlabel("$h$")
ylabel("error")
legend();
```
