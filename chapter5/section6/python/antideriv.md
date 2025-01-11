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
[**Demo %s**](#demo-int-antideriv)

The antiderivative of $e^x$ is, of course, itself. That makes evaluation of $\int_0^1 e^x\,dx$ by the Fundamental Theorem trivial.

```{code-cell}
exact = exp(1) - 1
```

```{index} ! Julia; quadgk
```

The module `scipy.integrate` has multiple functions that estimate the value of an integral numerically without finding the antiderivative first. As you can see here, it's often just as accurate.


```{code-cell}
from scipy.integrate import quad
Q, errest = quad(exp, 0, 1, epsabs=1e-13, epsrel=1e-13)
print(Q)

```

The numerical approach is also far more robust. For example, $e^{\,\sin x}$ has no useful antiderivative. But numerically, it's no more difficult.

```{code-cell}
Q, errest = quad(lambda x: exp(sin(x)), 0, 1, epsabs=1e-13, epsrel=1e-13)
print(Q)
```

When you look at the graphs of these functions, what's remarkable is that one of these areas is basic calculus while the other is almost impenetrable analytically. From a numerical standpoint, they are practically the same problem.

```{code-cell}
x = linspace(0, 1, 300)
subplot(1, 2, 1)
plot(x, exp(x))
ylim([0, 2.7]), title("exp(x)")
subplot(1, 2, 2)
plot(x, exp(sin(x)))
ylim([0, 2.7]), title("exp(sin(x))");
```
