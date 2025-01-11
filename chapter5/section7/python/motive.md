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
[**Demo %s**](#demo-adapt-motive)

This function gets increasingly oscillatory as $x$ increases.

```{code-cell}
f = lambda x: (x + 1) ** 2 * cos((2 * x + 1) / (x - 4.3))
x = linspace(0, 4, 600)
plot(x, f(x))
xlabel("$x$")
ylabel("$f(x)$");
```

Accordingly, the trapezoid rule is more accurate on the left half of this interval than on the right half.

```{code-cell}
n_ = 50 * 2 ** arange(4)
Tleft = zeros(4)
Tright = zeros(4)
for i, n in enumerate(n_):
    Tleft[i] = FNC.trapezoid(f, 0, 2, n)[0]
    Tright[i] = FNC.trapezoid(f, 2, 4, n)[0]
print("left half:", Tleft)
print("right half:", Tright)
```

```{code-cell}
from scipy.integrate import quad
left_val, err = quad(f, 0, 2, epsabs=1e-13, epsrel=1e-13)
right_val, err = quad(f, 2, 4, epsabs=1e-13, epsrel=1e-13)

print("    n     left error   right error")
for k in range(n_.size):
    print(f"  {n_[k]:4}    {Tleft[k]-left_val:8.3e}    {Tright[k]-right_val:8.3e}")
```

Both the picture and the numerical results suggest that more nodes should be used on the right half of the interval than on the left half.
