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
[**Demo %s**](#demo-int-extrap)

We estimate $\displaystyle\int_0^2 x^2 e^{-2x}\, dx$ using extrapolation. First we use `quadgk` to get an accurate value.

```{code-cell}
from scipy.integrate import quad
f = lambda x: x**2 * exp(-2 * x)
a = 0
b = 2
I, errest = quad(f, a, b, epsabs=1e-13, epsrel=1e-13)
print(f"Integral = {I:.14f}")
```

We start with the trapezoid formula on $n=N$ nodes.

```{code-cell}
N = 20    # the coarsest formula
n = N
h = (b - a) / n
t = h * arange(n + 1)
y = f(t)
```

We can now apply weights to get the estimate $T_f(N)$.

```{code-cell}
T = zeros(3)
T[0] = h * (sum(y[1:-1]) + y[0] / 2 + y[-1] / 2)
print(f"error (2nd order): {I - T[0]:.2e}")
```

Now we double to $n=2N$, but we only need to evaluate $f$ at every other interior node and apply {eq}`nc-doubling`.

```{code-cell}
n = 2 * n
h = h / 2
t = h * arange(n + 1)
T[1] = T[0] / 2 + h * sum(f(t[1:-1:2]))
print("error (2nd order):", I - T[:2])
```

As expected for a second-order estimate, the error went down by a factor of about 4. We can repeat the same code to double $n$ again.

```{code-cell}
n = 2 * n
h = h / 2
t = h * arange(n + 1)
T[2] = T[1] / 2 + h * sum(f(t[1:-1:2]))
print("error (2nd order):", I - T[:3])
```

Let us now do the first level of extrapolation to get results from Simpson's formula. We combine the elements `T[i]` and `T[i+1]` the same way for $i=1$ and $i=2$.

```{code-cell}
S = array([(4 * T[i + 1] - T[i]) / 3 for i in range(2)])
print("error (4th order):", I - S)
```

With the two Simpson values $S_f(N)$ and $S_f(2N)$ in hand, we can do one more level of extrapolation to get a sixth-order accurate result.

```{code-cell}
R = (16 * S[1] - S[0]) / 15
print("error (6th order):", I - R)
```

We can make a triangular table of the errors:

```{code-cell}
err = nan * ones((3, 3))
err[0, :] = I - T
err[1, 1:] = I - S
err[2, 2] = I - R
results = PrettyTable(["2nd order", "4th order", "6th order"])
results.add_rows(err.T)
print(results)
```

If we consider the computational time to be dominated by evaluations of $f$, then we have obtained a result with about twice as many accurate digits as the best trapezoid result, at virtually no extra cost.
