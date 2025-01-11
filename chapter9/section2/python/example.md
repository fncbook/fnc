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
[**Demo %s**](#demo-barycentric-example)

```{code-cell}
f = lambda x: sin(exp(2 * x))
x = linspace(0, 1, 500)
fig, ax = subplots()
ax.plot(x, f(x), label="function")
```

```{code-cell}
t = linspace(0, 1, 4)
y = f(t)
p = FNC.polyinterp(t, y)

ax.plot(x, p(x), label="interpolant")
ax.plot(t, y, "ko", label="nodes")
ax.legend()
ax.set_title("Interpolation on 4 nodes")
fig
```

The curves must intersect at the interpolation nodes. For $n=6$ the interpolant is noticeably better.

```{code-cell}
plot(x, f(x), label="function")
t = linspace(0, 1, 7)
y = f(t)
p = FNC.polyinterp(t, y)
plot(x, p(x), label="interpolant")
plot(t, y, "ko", label="nodes")
legend(),  title("Interpolation on 7 nodes");
```
