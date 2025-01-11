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
[**Demo %s**](#demo-finitediffs-fd1)

If $f(x)=e^{\,\sin(x)}$, then $f'(0)=1$.

```{code-cell}
f = lambda x: exp(sin(x))
```

Here are the first two centered differences from {numref}`table-FDcenter`.

```{code-cell}
h = 0.05
CD2 = (-f(-h) + f(h)) / (2*h)
CD4 = (f(-2*h) - 8*f(-h) + 8*f(h) - f(2*h)) / (12*h)
print(f"CD2 is {CD2:.9f} and CD4 is {CD4:.9f}")
```

Here are the first two forward differences from {numref}`table-FDforward`.

```{code-cell}
FD1 = (-f(0) + f(h)) / h
FD2 = (-3*f(0) + 4*f(h) - f(2*h)) / (2*h)
print(f"FD1 is {FD1:.9f} and FD2 is {FD2:.9f}")
```

Finally, here are the backward differences that come from reverse-negating the forward differences.

```{code-cell}
BD1 = (-f(-h) + f(0)) / h
BD2 = (f(-2*h) - 4*f(-h) + 3*f(0)) / (2*h)
print(f"BD1 is {BD1:.9f} and BD2 is {BD2:.9f}")
```
