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
[**Demo %s**](#demo-newton-usage)

Suppose we want to evaluate the inverse of the function $h(x)=e^x-x$. This means solving $y=h(x)$, or $h(x)-y=0$, for $x$ when $y$ is given. That equation has no solution in terms of elementary functions. If a value of $y$ is given numerically, though, we simply have a rootfinding problem for $f(x)=e^x-x-y$.
```{tip}
:class: dropdown
The `enumerate` function produces a pair of values for each iteration: a positional index and the corresponding contents.
```

```{index} ! Python; enumerate
```

```{code-cell}
h = lambda x: exp(x) - x
dh_dx = lambda x: exp(x) - 1
y_ = linspace(h(0), h(2), 200)
x_ = zeros(y_.shape)
for (i, y) in enumerate(y_):
    f = lambda x: h(x) - y
    df_dx = lambda x: dh_dx(x)
    x = FNC.newton(f, df_dx, y)
    x_[i] = x[-1]

plot(x_, y_, label="$y=h(x)$")
plot(y_, x_, label="$y=h^{-1}(x)$")
plot([0, max(y_)], [0, max(y_)], 'k--', label="")
title("Function and its inverse")
xlabel("x"), ylabel("y"), axis("equal")
ax.grid()
legend();
```
