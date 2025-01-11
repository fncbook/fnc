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
[**Demo %s**](#demo-finitediffs-fd2)

If $f(x)=e^{\,\sin(x)}$, then $f''(0)=1$.

```{code-cell}
f = lambda x: exp(sin(x))
```

Here is a centered estimate given by {eq}`centerFD22`.

```{code-cell}
h = 0.05
CD2 = (f(-h) - 2*f(0) + f(h)) / h**2
print(f"CD2 is {CD2:.9f}")
```

For the same $h$, here are forward estimates given by {eq}`forwardFD21` and {eq}`forwardFD22`.

```{code-cell}
FD1 = (f(0) - 2*f(h) + f(2*h)) / h**2
FD2 = (2*f(0) - 5*f(h) + 4*f(2*h) - f(3*h)) / h**2
print(f"FD1 is {FD1:.9f} and FD2 is {FD2:.9f}")
```

Finally, here are the backward estimates that come from reversing {eq}`forwardFD21` and {eq}`forwardFD22`.

```{code-cell}
BD1 = (f(-2*h) - 2*f(-h) + f(0)) / h**2
BD2 = (-f(-3*h) + 4*f(-2*h) - 5*f(-h) + 2*f(0)) / h**2
print(f"BD1 is {BD1:.9f} and BD2 is {BD2:.9f}")
```
