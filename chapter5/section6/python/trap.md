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
[**Demo %s**](#demo-int-trap)

We will approximate the integral of the function $f(x)=e^{\sin 7x}$ over the interval $[0,2]$.

```{code-cell}
f = lambda x: exp(sin(7 * x))
a, b = 0, 2
```

In lieu of the exact value, we will use the `quad` function to find an accurate result.

```{code-cell}
from scipy.integrate import quad
I, errest = quad(f, a, b, epsabs=1e-13, epsrel=1e-13)
print(f"Integral = {I:.14f}")
```

Here is the trapezoid result at $n=40$, and its error.

```{code-cell}
T, t, y = FNC.trapezoid(f, a, b, 40)
print(f"Trapezoid estimate is {T:.14f} with error {I - T:.2e}")
```

In order to check the order of accuracy, we increase $n$ by orders of magnitude and observe how the error decreases.

```{code-cell}
n_ = 40 * 2 ** arange(6)
err = zeros(size(n_))
print("     n     error")
for k, n in enumerate(n_):
    T, t, y = FNC.trapezoid(f, a, b, n)
    err[k] = I - T
    print(f"{n:6d}   {err[k]:8.3e} ")
```

Each increase by a factor of 10 in $n$ cuts the error by a factor of about 100, which is consistent with second-order convergence. Another check is that a log-log graph should give a line of slope $-2$ as $n\to\infty$.

```{code-cell}
loglog(n_, abs(err), "-o", label="results")
loglog(n_, 3e-3 * (n_ / n_[0]) ** (-2), "--", label="2nd order")
gca().invert_xaxis()
xlabel("$n$")
ylabel("error")
legend()
title("Convergence of trapezoidal integration");
